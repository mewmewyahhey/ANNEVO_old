import unittest
from unittest.mock import patch

import numpy as np

from src.gene_decoding import get_gene_region, process_gene_segment


class _FakeSeqRecord:
    def __init__(self, seq):
        self.seq = seq


class GeneDecodingWorkerTest(unittest.TestCase):
    @patch("src.gene_decoding.h5py.File", side_effect=AssertionError("worker should not open H5"))
    @patch("src.gene_decoding.decode_gene_structure", return_value=[("gene",)])
    def test_process_gene_segment_uses_supplied_payload(self, mock_decode, mock_h5):
        prediction_slice = np.ones((3, 5), dtype=np.float32)
        region = (0, 3, "chr1", 1, prediction_slice, "ATG")

        result = process_gene_segment(region, 60, 0.5, 1, 0)

        self.assertEqual(result, ([("gene",)], "chr1", 1))
        mock_h5.assert_not_called()
        mock_decode.assert_called_once()

    @patch("src.gene_decoding.detect_gene_location", return_value=[(1, 4)])
    def test_get_gene_region_supplies_in_memory_payload(self, mock_detect):
        predictions_forward = np.arange(30, dtype=np.float32).reshape(6, 5)
        predictions_reverse = np.arange(30, 60, dtype=np.float32).reshape(6, 5)
        genome_predictions = {"chr1": [predictions_forward, predictions_reverse]}
        genome_seq = {"chr1": _FakeSeqRecord("AACCGG")}

        _, regions = get_gene_region(genome_predictions, genome_seq, 0.1, 0.5)

        self.assertEqual(len(regions), 2)
        forward_region = regions[0]
        reverse_region = regions[1]
        np.testing.assert_array_equal(forward_region[4], predictions_forward[1:4])
        self.assertEqual(forward_region[5], "ACC")
        np.testing.assert_array_equal(reverse_region[4], predictions_reverse[1:4])
        self.assertEqual(reverse_region[5], "CGG")


if __name__ == "__main__":
    unittest.main()
