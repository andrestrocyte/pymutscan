import unittest

from pymutscan.collapse import group_similar_sequences, hamming_distance


class CollapseTests(unittest.TestCase):
    def test_mutscan_reference_examples(self):
        seqs = ["AACGTAGCA", "ACCGTAGCA", "AACGGAGCA", "ATCGGAGCA", "TGAGGCATA"]
        scores = [5, 1, 3, 1, 8]
        expected = [seqs[i] for i in [0, 0, 0, 3, 4]]
        mapping = group_similar_sequences(seqs, scores, collapse_max_dist=1)
        self.assertEqual([mapping[s] for s in seqs], expected)

        expected_ratio = [seqs[i] for i in [0, 0, 2, 2, 4]]
        mapping = group_similar_sequences(
            seqs, scores, collapse_max_dist=1, collapse_min_ratio=2
        )
        self.assertEqual([mapping[s] for s in seqs], expected_ratio)

    def test_tie_break_is_lexicographic(self):
        seqs = ["AAAT", "AAAA"]
        mapping = group_similar_sequences(seqs, [2, 2], collapse_max_dist=1)
        self.assertEqual(mapping, {"AAAT": "AAAA", "AAAA": "AAAA"})

    def test_hamming_requires_equal_length(self):
        with self.assertRaises(ValueError):
            hamming_distance("A", "AA")


if __name__ == "__main__":
    unittest.main()

