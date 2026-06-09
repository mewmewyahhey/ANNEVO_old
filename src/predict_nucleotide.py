from src.utils import model_construction, model_load_weights
from src.chunking_utils import contig_window_count, should_flush_prediction_chunk
from Bio import SeqIO
import gc
import h5py
import numpy as np
import time
import torch
import torch.nn.functional as F
from tqdm import tqdm
from torch.utils.data import DataLoader
from torch.utils.data import Dataset


WRITE_BUFFER_BATCHES = 4


def reverse_complement(dna_sequence):
    complement_map = str.maketrans('ATGCRMYWKBSHDVNXatgcrmywkbshdvnx', 'TACGRMYWKBSHDVNXtacgrmywkbshdvnx')
    return dna_sequence.translate(complement_map)[::-1]


class GenomeDataset(Dataset):
    def __init__(self, genome_data, reverse=False, index_order=None):
        self.data = genome_data
        self.reverse = reverse
        self.index_order = index_order

    def __len__(self):
        if self.index_order is None:
            return len(self.data)
        return len(self.index_order)

    def __getitem__(self, idx):
        data_idx = idx if self.index_order is None else self.index_order[idx]
        window_seq = self.data[data_idx]
        if self.reverse:
            window_seq = reverse_complement(window_seq)
        one_hot_seq = sequence_encode(window_seq)
        return torch.tensor(one_hot_seq, dtype=torch.float)


def sequence_encode(seq):
    mapping = {'A': [1, 0, 0, 0],
               'C': [0, 1, 0, 0],
               'G': [0, 0, 1, 0],
               'T': [0, 0, 0, 1],
               'N': [0.25, 0.25, 0.25, 0.25],
               'M': [0.25, 0.25, 0.25, 0.25],
               'W': [0.25, 0.25, 0.25, 0.25],
               'R': [0.25, 0.25, 0.25, 0.25],
               'Y': [0.25, 0.25, 0.25, 0.25],
               'K': [0.25, 0.25, 0.25, 0.25],
               'B': [0.25, 0.25, 0.25, 0.25],
               'S': [0.25, 0.25, 0.25, 0.25],
               'D': [0.25, 0.25, 0.25, 0.25],
               'H': [0.25, 0.25, 0.25, 0.25],
               'V': [0.25, 0.25, 0.25, 0.25],
               'X': [0, 0, 0, 0]}
    return [mapping[s] for s in seq]


def build_window_sequence(sequence_forward, start, window_size, flank_length):
    length = len(sequence_forward)
    end = start + window_size
    if start - flank_length < 0:
        if end + flank_length <= length:
            pad_before = 'X' * (flank_length - start)
            return pad_before + sequence_forward[0:end + flank_length]
        pad_before = 'X' * (flank_length - start)
        pad_after = 'X' * (end + flank_length - length)
        return pad_before + sequence_forward[0:length] + pad_after
    if end + flank_length > length:
        pad_after = 'X' * (end + flank_length - length)
        return sequence_forward[start - flank_length:length] + pad_after
    return sequence_forward[start - flank_length:end + flank_length]


def append_sequence_windows(sequence_forward, window_size, flank_length, windows_forward, reverse_window_indices):
    sequence_window_start = len(windows_forward)
    for start in range(0, len(sequence_forward), window_size):
        windows_forward.append(build_window_sequence(sequence_forward, start, window_size, flank_length))
    sequence_window_end = len(windows_forward)
    # Reverse-strand prediction consumes each sequence's windows in reverse order.
    reverse_window_indices.extend(range(sequence_window_end - 1, sequence_window_start - 1, -1))
    return sequence_window_end - sequence_window_start


def create_sequence_layout(seq_id_chunk, seq_length_chunk, offset, window_size):
    seq_layout = []
    for i, seq_id in enumerate(seq_id_chunk):
        length = seq_length_chunk[i]
        global_start = offset[i] * window_size
        global_end = offset[i + 1] * window_size
        seq_layout.append({
            "seq_id": seq_id,
            "forward_data_start": global_start,
            "forward_data_end": global_start + length,
            "reverse_data_start": global_end - length,
            "reverse_data_end": global_end,
        })
    return seq_layout


def create_prediction_datasets(h5file, seq_layout, num_classes):
    for seq_info in seq_layout:
        seq_id = seq_info["seq_id"]
        if seq_id in h5file:
            raise ValueError(f"Duplicate sequence id '{seq_id}' found while writing prediction output.")
        group = h5file.create_group(seq_id)
        forward_length = seq_info["forward_data_end"] - seq_info["forward_data_start"]
        reverse_length = seq_info["reverse_data_end"] - seq_info["reverse_data_start"]
        group.create_dataset("predictions_forward", shape=(forward_length, num_classes), dtype='float16')
        group.create_dataset("predictions_reverse", shape=(reverse_length, num_classes), dtype='float16')


def write_prediction_block(h5file, seq_layout, dataset_name, data_start_key, data_end_key, block_predictions, block_start, seq_index):
    block_end = block_start + len(block_predictions)

    while seq_index < len(seq_layout) and seq_layout[seq_index][data_end_key] <= block_start:
        seq_index += 1

    # Scatter the flattened prediction block back into per-sequence HDF5 datasets.
    current_index = seq_index
    while current_index < len(seq_layout):
        seq_info = seq_layout[current_index]
        data_start = seq_info[data_start_key]
        if data_start >= block_end:
            break

        data_end = seq_info[data_end_key]
        overlap_start = max(block_start, data_start)
        overlap_end = min(block_end, data_end)
        if overlap_start < overlap_end:
            block_slice = slice(overlap_start - block_start, overlap_end - block_start)
            dataset_slice = slice(overlap_start - data_start, overlap_end - data_start)
            h5file[seq_info["seq_id"]][dataset_name][dataset_slice] = block_predictions[block_slice]

        if data_end <= block_end:
            current_index += 1
        else:
            break

    while seq_index < len(seq_layout) and seq_layout[seq_index][data_end_key] <= block_end:
        seq_index += 1

    return seq_index


def flush_prediction_blocks(h5file, seq_layout, dataset_name, data_start_key, data_end_key, pending_blocks, block_start, seq_index):
    if not pending_blocks:
        return 0.0, block_start, seq_index

    start_time = time.time()
    if len(pending_blocks) == 1:
        block_predictions = pending_blocks[0]
    else:
        block_predictions = np.concatenate(pending_blocks, axis=0)
    seq_index = write_prediction_block(
        h5file,
        seq_layout,
        dataset_name,
        data_start_key,
        data_end_key,
        block_predictions,
        block_start,
        seq_index,
    )
    runtime = time.time() - start_time
    block_start += len(block_predictions)
    pending_blocks.clear()
    return runtime, block_start, seq_index


def stream_predictions(model, windows_forward, reverse_window_indices, device, num_classes, batch_size, num_workers,
                       seq_layout, h5file, window_size, reverse=False):
    if reverse:
        data = GenomeDataset(windows_forward, reverse=True, index_order=reverse_window_indices)
        dataset_name = "predictions_reverse"
        data_start_key = "reverse_data_start"
        data_end_key = "reverse_data_end"
    else:
        data = GenomeDataset(windows_forward)
        dataset_name = "predictions_forward"
        data_start_key = "forward_data_start"
        data_end_key = "forward_data_end"

    dataloader = DataLoader(data, batch_size=batch_size, shuffle=False, num_workers=num_workers)
    write_buffer_rows = max(window_size * max(batch_size, 1) * WRITE_BUFFER_BATCHES, window_size)
    pending_blocks = []
    pending_rows = 0
    block_start = 0
    seq_index = 0
    file_saving_time = 0.0

    with torch.no_grad():
        for data in tqdm(dataloader):
            seqs = data.to(device)
            outputs, _, _ = model(seqs)
            outputs = outputs.reshape(-1, num_classes)
            outputs = F.softmax(outputs, dim=-1).cpu().numpy().astype('float16')
            pending_blocks.append(outputs)
            pending_rows += outputs.shape[0]
            del seqs

            if pending_rows >= write_buffer_rows:
                runtime, block_start, seq_index = flush_prediction_blocks(
                    h5file,
                    seq_layout,
                    dataset_name,
                    data_start_key,
                    data_end_key,
                    pending_blocks,
                    block_start,
                    seq_index,
                )
                file_saving_time += runtime
                pending_rows = 0

        runtime, block_start, seq_index = flush_prediction_blocks(
            h5file,
            seq_layout,
            dataset_name,
            data_start_key,
            data_end_key,
            pending_blocks,
            block_start,
            seq_index,
        )
        file_saving_time += runtime

    return file_saving_time


def predict_chunk(model, windows_forward, reverse_window_indices, device, num_classes, batch_size, num_workers,
                  seq_id_chunk, seq_length_chunk, offset, window_size, h5file):
    seq_layout = create_sequence_layout(seq_id_chunk, seq_length_chunk, offset, window_size)
    create_prediction_datasets(h5file, seq_layout, num_classes)
    file_saving_time = stream_predictions(
        model,
        windows_forward,
        reverse_window_indices,
        device,
        num_classes,
        batch_size,
        num_workers,
        seq_layout,
        h5file,
        window_size,
        reverse=False,
    )
    file_saving_time += stream_predictions(
        model,
        windows_forward,
        reverse_window_indices,
        device,
        num_classes,
        batch_size,
        num_workers,
        seq_layout,
        h5file,
        window_size,
        reverse=True,
    )
    return file_saving_time


def nucleotide_prediction(genome, model_path, genome_size_threshold, num_workers, prediction_path, batch_size, max_windows_per_chunk, window_size, flank_length, channels, dim_feedforward,
                          num_encoder_layers, num_heads, num_blocks, num_branches, num_classes):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = model_construction(device, window_size, flank_length, channels, dim_feedforward, num_encoder_layers, num_heads, num_blocks, num_branches, num_classes, top_k=2)
    model = model_load_weights(model_path, model, device)
    model.eval()
    chunk_num = 1

    file_saving_time = 0

    windows_forward = []
    reverse_window_indices = []
    offset = [0]
    count = 0
    cumulative_size = 0
    seq_id_chunk = []
    seq_length_chunk = []
    with h5py.File(prediction_path, "w") as h5file:
        with open(genome) as fna:
            for chromosome_seq_record in SeqIO.parse(fna, "fasta"):
                chromosome = chromosome_seq_record.id
                sequence_forward = str(chromosome_seq_record.seq).upper()
                length = len(sequence_forward)
                chromosome_windows = contig_window_count(length, window_size)

                if should_flush_prediction_chunk(
                    chunk_has_content=bool(seq_id_chunk),
                    current_windows=count,
                    next_contig_windows=chromosome_windows,
                    max_windows_per_chunk=max_windows_per_chunk,
                ):
                    print(f'---------------------------------------Prediction on chunk {chunk_num}---------------------------------------')
                    chunk_num += 1
                    runtime = predict_chunk(
                        model,
                        windows_forward,
                        reverse_window_indices,
                        device,
                        num_classes,
                        batch_size,
                        num_workers,
                        seq_id_chunk,
                        seq_length_chunk,
                        offset,
                        window_size,
                        h5file,
                    )
                    file_saving_time += runtime

                    windows_forward = []
                    reverse_window_indices = []
                    offset = [0]
                    count = 0
                    cumulative_size = 0
                    seq_id_chunk = []
                    seq_length_chunk = []
                    gc.collect()
                    if device.type == "cuda":
                        torch.cuda.empty_cache()

                cumulative_size += length
                seq_id_chunk.append(chromosome)
                seq_length_chunk.append(length)
                count += append_sequence_windows(sequence_forward, window_size, flank_length, windows_forward, reverse_window_indices)
                offset.append(count)

                if cumulative_size > genome_size_threshold:
                    print(f'---------------------------------------Prediction on chunk {chunk_num}---------------------------------------')
                    chunk_num += 1
                    runtime = predict_chunk(
                        model,
                        windows_forward,
                        reverse_window_indices,
                        device,
                        num_classes,
                        batch_size,
                        num_workers,
                        seq_id_chunk,
                        seq_length_chunk,
                        offset,
                        window_size,
                        h5file,
                    )
                    file_saving_time += runtime

                    windows_forward = []
                    reverse_window_indices = []
                    offset = [0]
                    count = 0
                    cumulative_size = 0
                    seq_id_chunk = []
                    seq_length_chunk = []
                    gc.collect()
                    if device.type == "cuda":
                        torch.cuda.empty_cache()

        if seq_id_chunk:
            print(f'---------------------------------------Prediction on chunk {chunk_num}---------------------------------------')
            chunk_num += 1
            runtime = predict_chunk(
                model,
                windows_forward,
                reverse_window_indices,
                device,
                num_classes,
                batch_size,
                num_workers,
                seq_id_chunk,
                seq_length_chunk,
                offset,
                window_size,
                h5file,
            )
            file_saving_time += runtime

    print(f"file saving cost {file_saving_time:.1f} seconds")

    del windows_forward
    del reverse_window_indices
    del seq_id_chunk
    del seq_length_chunk
    del model
    torch.cuda.empty_cache()
    gc.collect()
