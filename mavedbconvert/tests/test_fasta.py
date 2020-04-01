import os
import unittest

from mavedbconvert.tests import ProgramTestCase

from mavedbconvert.fasta import parse_fasta, split_fasta_path


class TestFastaPath(ProgramTestCase):
    def test_infers_bzip(self):
        head, base, ext, compression = split_fasta_path(
            os.path.join(self.data_dir, "fasta", "wt.fasta.bz2")
        )
        self.assertEqual(ext, ".fasta")
        self.assertEqual(compression, "bz2")

    def test_infers_gzip(self):
        head, base, ext, compression = split_fasta_path(
            os.path.join(self.data_dir, "fasta", "wt.fasta.gz")
        )
        self.assertEqual(ext, ".fasta")
        self.assertEqual(compression, "gz")

    def test_infers_uncompressed(self):
        head, base, ext, compression = split_fasta_path(
            os.path.join(self.data_dir, "fasta", "wt.fasta")
        )
        self.assertEqual(ext, ".fasta")
        self.assertEqual(compression, None)

        head, base, ext, compression = split_fasta_path(
            os.path.join(self.data_dir, "fasta", "lower.fa")
        )
        self.assertEqual(ext, ".fa")
        self.assertEqual(compression, None)

    def test_ioerror_invalid_ext(self):
        with self.assertRaises(IOError):
            split_fasta_path(os.path.join(self.data_dir, "enrich", "enrich.tsv"))


class TestFastaReader(ProgramTestCase):
    def test_can_read_first_sequence(self):
        sequence = parse_fasta(os.path.join(self.data_dir, "fasta", "wt.fasta"))
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTAA"
        )
        self.assertEqual(sequence, expected)

    def test_converts_to_uppercase(self):
        sequence = parse_fasta(os.path.join(self.data_dir, "fasta", "lower.fa"))
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTA"
        )
        self.assertEqual(sequence, expected)

    def test_error_more_than_one_sequence(self):
        with self.assertRaises(ValueError):
            parse_fasta(os.path.join(self.data_dir, "fasta", "two.fasta"))

    def test_error_invalid_chars_in_sequence(self):
        with self.assertRaises(ValueError):
            parse_fasta(os.path.join(self.data_dir, "fasta", "invalid_chars.fasta"))

    def test_ignores_blank_lines(self):
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTAA"
        )
        seq = parse_fasta(os.path.join(self.data_dir, "fasta", "spaces.fasta"))
        self.assertEqual(seq, expected)

    def test_error_missing_gt_on_first_line(self):
        with self.assertRaises(IOError):
            parse_fasta(os.path.join(self.data_dir, "fasta", "bad_format.fasta"))

    def test_can_open_with_gzip(self):
        sequence = parse_fasta(os.path.join(self.data_dir, "fasta", "wt.fasta.gz"))
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTAA"
        )
        self.assertEqual(sequence, expected)

    def test_can_open_with_bzip(self):
        sequence = parse_fasta(os.path.join(self.data_dir, "fasta", "wt.fasta.bz2"))
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTAA"
        )
        self.assertEqual(sequence, expected)


if __name__ == "__main__":
    unittest.main()
