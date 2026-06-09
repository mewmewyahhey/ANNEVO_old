import unittest

from src.chunking_utils import contig_window_count, should_flush_prediction_chunk


class ChunkingUtilsTest(unittest.TestCase):
    def test_contig_window_count_rounds_up(self):
        self.assertEqual(contig_window_count(1, 30720), 1)
        self.assertEqual(contig_window_count(30720, 30720), 1)
        self.assertEqual(contig_window_count(30721, 30720), 2)

    def test_window_cap_flushes_before_adding_next_contig(self):
        self.assertTrue(
            should_flush_prediction_chunk(
                chunk_has_content=True,
                current_windows=8000,
                next_contig_windows=300,
                max_windows_per_chunk=8192,
            )
        )

    def test_single_oversized_contig_is_allowed_to_form_its_own_chunk(self):
        self.assertFalse(
            should_flush_prediction_chunk(
                chunk_has_content=False,
                current_windows=0,
                next_contig_windows=9000,
                max_windows_per_chunk=8192,
            )
        )

    def test_zero_or_negative_window_cap_disables_limit(self):
        self.assertFalse(
            should_flush_prediction_chunk(
                chunk_has_content=True,
                current_windows=8000,
                next_contig_windows=300,
                max_windows_per_chunk=0,
            )
        )


if __name__ == '__main__':
    unittest.main()
