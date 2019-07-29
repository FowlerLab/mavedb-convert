import os
from unittest import TestCase

from ..fasta import parse_fasta, split_fasta_path

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(BASE_DIR + "/data/")


class TestSplitFasta(TestCase):
    def test_infers_bzip(self):
        head, base, ext, compression = split_fasta_path(
            os.path.join(DATA_DIR, "wt.fasta.bz2")
        )
        self.assertEqual(ext, ".fasta")
        self.assertEqual(compression, "bz2")

    def test_infers_gzip(self):
        head, base, ext, compression = split_fasta_path(
            os.path.join(DATA_DIR, "wt.fasta.gz")
        )
        self.assertEqual(ext, ".fasta")
        self.assertEqual(compression, "gz")

    def test_infers_none(self):
        head, base, ext, compression = split_fasta_path(
            os.path.join(DATA_DIR, "wt.fasta")
        )
        self.assertEqual(ext, ".fasta")
        self.assertEqual(compression, None)

        head, base, ext, compression = split_fasta_path(
            os.path.join(DATA_DIR, "lower.fa")
        )
        self.assertEqual(ext, ".fa")
        self.assertEqual(compression, None)

    def test_ioerror_invalid_ext(self):
        with self.assertRaises(IOError):
            split_fasta_path(os.path.join(DATA_DIR, "enrich1.tsv"))


class TestFastaReader(TestCase):
    def test_can_read_first_sequence(self):
        sequence = parse_fasta(os.path.join(DATA_DIR, "wt.fasta"))
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTAA"
        )
        self.assertEqual(sequence, expected)

    def test_converts_to_uppercase(self):
        sequence = parse_fasta(os.path.join(DATA_DIR, "lower.fa"))
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTA"
        )
        self.assertEqual(sequence, expected)

    def test_error_more_than_one_sequence(self):
        with self.assertRaises(ValueError):
            parse_fasta(os.path.join(DATA_DIR, "two.fasta"))

    def test_error_invalid_chars_in_sequence(self):
        with self.assertRaises(ValueError):
            parse_fasta(os.path.join(DATA_DIR, "invalid_chars.fasta"))

    def test_ignores_blank_lines(self):
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTAA"
        )
        seq = parse_fasta(os.path.join(DATA_DIR, "spaces.fasta"))
        self.assertEqual(seq, expected)

    def test_error_missing_gt_on_first_line(self):
        with self.assertRaises(IOError):
            parse_fasta(os.path.join(DATA_DIR, "bad_format.fasta"))

    def test_can_open_with_gzip(self):
        sequence = parse_fasta(os.path.join(DATA_DIR, "wt.fasta.gz"))
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTAA"
        )
        self.assertEqual(sequence, expected)

    def test_can_open_with_bzip(self):
        sequence = parse_fasta(os.path.join(DATA_DIR, "wt.fasta.bz2"))
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTAA"
        )
        self.assertEqual(sequence, expected)
