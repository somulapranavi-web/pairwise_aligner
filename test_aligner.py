#!/usr/bin/env python3
"""
Unit tests for the Needleman-Wunsch Sequence Aligner.
"""

import unittest
from aligner import needleman_wunsch

class TestNeedlemanWunsch(unittest.TestCase):
    
    def test_identical_single_character(self):
        """Test alignment of identical single-character sequences."""
        score, aligned1, aligned2, _ = needleman_wunsch("A", "A", match=1, mismatch=-1, gap=-2)
        self.assertEqual(score, 1)
        self.assertEqual(aligned1, "A")
        self.assertEqual(aligned2, "A")

    def test_identical_sequences(self):
        """Test alignment of longer identical sequences."""
        score, aligned1, aligned2, _ = needleman_wunsch("GCAT", "GCAT", match=1, mismatch=-1, gap=-2)
        self.assertEqual(score, 4)
        self.assertEqual(aligned1, "GCAT")
        self.assertEqual(aligned2, "GCAT")

    def test_mismatch_only(self):
        """Test alignment of sequences that differ by mismatches (no gaps needed)."""
        # HELLO vs CELLO
        # H-C (mismatch: -1) + ELLO-ELLO (matches: 4) = 3
        score, aligned1, aligned2, _ = needleman_wunsch("HELLO", "CELLO", match=1, mismatch=-1, gap=-2)
        self.assertEqual(score, 3)
        self.assertEqual(aligned1, "HELLO")
        self.assertEqual(aligned2, "CELLO")

    def test_gap_handling(self):
        """Test that gaps are placed correctly when one sequence is longer."""
        # COLD vs FOLDED
        # C-F (mismatch: -1), O-O (match: 1), L-L (match: 1), D-D (match: 1), --E (gap: -2), --D (gap: -2)
        # Score: -1 + 1 + 1 + 1 - 2 - 2 = -2
        score, aligned1, aligned2, _ = needleman_wunsch("COLD", "FOLDED", match=1, mismatch=-1, gap=-2)
        self.assertEqual(score, -2)
        self.assertEqual(aligned1, "COL--D")
        self.assertEqual(aligned2, "FOLDED")

    def test_empty_sequences(self):
        """Test handling of empty sequence inputs."""
        # Both empty
        score, aligned1, aligned2, _ = needleman_wunsch("", "", match=1, mismatch=-1, gap=-2)
        self.assertEqual(score, 0)
        self.assertEqual(aligned1, "")
        self.assertEqual(aligned2, "")
        
        # One empty
        score, aligned1, aligned2, _ = needleman_wunsch("A", "", match=1, mismatch=-1, gap=-2)
        self.assertEqual(score, -2)
        self.assertEqual(aligned1, "A")
        self.assertEqual(aligned2, "-")

    def test_dna_sequences_with_gaps(self):
        """Test a standard DNA sequence alignment with gaps."""
        # G A T T A C A
        # G C A T - C A  (with gap at index 4)
        # G-G (1), A-C (-1), T-A (-1), T-T (1), A-- (-2), C-C (1), A-A (1)
        # Total score: 1 - 1 - 1 + 1 - 2 + 1 + 1 = 0
        score, aligned1, aligned2, _ = needleman_wunsch("GATTACA", "GCATCA", match=1, mismatch=-1, gap=-2)
        self.assertEqual(score, 0)
        self.assertEqual(aligned1, "GATTACA")
        self.assertEqual(aligned2, "GCAT-CA")

class TestSmithWaterman(unittest.TestCase):
    
    def test_local_substring(self):
        """Test that it finds a matching substring in different contexts."""
        from aligner import smith_waterman
        score, aligned1, aligned2, _ = smith_waterman(
            "CCCCCGATTACACCCCC", 
            "GGGGGGATTACAGGGGGGG", 
            match=1, mismatch=-1, gap=-2
        )
        self.assertEqual(score, 7)
        self.assertEqual(aligned1, "GATTACA")
        self.assertEqual(aligned2, "GATTACA")

    def test_no_alignment(self):
        """Test that completely different sequences yield an empty local alignment and 0 score."""
        from aligner import smith_waterman
        score, aligned1, aligned2, _ = smith_waterman(
            "AAAAA", 
            "CCCCC", 
            match=1, mismatch=-1, gap=-2
        )
        self.assertEqual(score, 0)
        self.assertEqual(aligned1, "")
        self.assertEqual(aligned2, "")

if __name__ == "__main__":
    unittest.main()
