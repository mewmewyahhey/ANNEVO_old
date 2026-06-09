DEFAULT_MAX_WINDOWS_PER_CHUNK = 8192


def contig_window_count(sequence_length, window_size):
    if sequence_length <= 0:
        return 0
    return (sequence_length + window_size - 1) // window_size


def should_flush_prediction_chunk(chunk_has_content, current_windows, next_contig_windows, max_windows_per_chunk):
    if not chunk_has_content:
        return False
    if max_windows_per_chunk is None or max_windows_per_chunk <= 0:
        return False
    return current_windows + next_contig_windows > max_windows_per_chunk
