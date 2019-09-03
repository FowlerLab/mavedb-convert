import os
import mock
from unittest import TestCase

from .. import parsers, exceptions, constants

from . import ProgramTestCase


TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


class TestParseBoolean(TestCase):
    def test_true_if_str_of_true(self):
        self.assertTrue(parsers.parse_boolean(True))
        self.assertTrue(parsers.parse_boolean("True"))

    def test_false_if_not_repr_of_true(self):
        self.assertFalse(parsers.parse_boolean(None))
        self.assertFalse(parsers.parse_boolean("none"))
        self.assertFalse(parsers.parse_boolean(""))
        self.assertFalse(parsers.parse_boolean(False))


class TestParseNumeric(TestCase):
    def test_converts_to_dtype(self):
        self.assertIsInstance(
            parsers.parse_numeric("1", name="int", dtype=float), float
        )

    def test_value_error_cannot_cast_to_dtype(self):
        with self.assertRaises(ValueError):
            parsers.parse_numeric(None, name="value", dtype=int)
        with self.assertRaises(ValueError):
            parsers.parse_numeric("a", name="value", dtype=int)


class TestParseString(TestCase):
    def test_returns_none_if_falsey(self):
        self.assertIsNone(parsers.parse_string(None))
        self.assertIsNone(parsers.parse_string(" "))
        self.assertIsNone(parsers.parse_string(""))

    def test_returns_string_stripped_of_ws(self):
        self.assertEqual(parsers.parse_string(" aaa "), "aaa")


class TestParseSrc(TestCase):
    @mock.patch(
        "mavedbconvert.parsers.parse_string",
        return_value=os.path.join(TEST_DATA_DIR, "enrich2", "enrich2.tsv"),
    )
    def test_calls_parse_string(self, patch):
        parsers.parse_src(os.path.join(TEST_DATA_DIR, "enrich2", "enrich2.tsv"))
        patch.assert_called()

    def test_ok_file_exists(self):
        path = os.path.join(TEST_DATA_DIR, "enrich2", "enrich2.tsv")
        self.assertEqual(path, parsers.parse_src(path))

    def test_error_no_value(self):
        for v in (None, "None", ""):
            with self.assertRaises(ValueError):
                parsers.parse_src(v)

    def test_error_file_not_found(self):
        path = os.path.join(TEST_DATA_DIR, "enrich2", "missing_file.tsv")
        with self.assertRaises(FileNotFoundError):
            parsers.parse_src(path)

    def test_error_file_is_a_dir(self):
        with self.assertRaises(IsADirectoryError):
            parsers.parse_src(os.path.join(TEST_DATA_DIR))


class TestParseDst(ProgramTestCase):
    @mock.patch(
        "mavedbconvert.parsers.parse_string", return_value=os.path.join(TEST_DATA_DIR)
    )
    def test_calls_parse_string(self, patch):
        parsers.parse_dst(os.path.join(TEST_DATA_DIR))
        patch.assert_called()

    def test_ok_dst_exists(self):
        path = os.path.join(os.path.join(TEST_DATA_DIR))
        self.assertEqual(path, parsers.parse_dst(path))

    def test_returns_none_no_value(self):
        for v in (None, "None", ""):
            self.assertIsNone(parsers.parse_dst(v))

    def test_dst_path_is_normalised(self):
        path = TEST_DATA_DIR + "//fasta"
        self.assertEqual(parsers.parse_dst(path), os.path.join(TEST_DATA_DIR, "fasta"))

    def test_makes_dst_directory_tree(self):
        path = os.path.join(TEST_DATA_DIR, "subdir")
        parsers.parse_dst(path)
        self.assertTrue(os.path.isdir(path))
        self.bin.append(path)


class TestParseProgram(TestCase):
    def test_ok_supported_program(self):
        for p in ("enrich2", "enrich", "empiric"):
            parsers.parse_program(p)

    def test_error_unsupported_program(self):
        with self.assertRaises(ValueError):
            parsers.parse_program("aaa")

    def test_sets_correct_program_from_dict(self):
        program = {"enrich": True, "empiric": False, "enrich2": False}
        self.assertEqual(parsers.parse_program(program), "enrich")

        program = {"enrich": False, "empiric": False, "enrich2": True}
        self.assertEqual(parsers.parse_program(program), "enrich2")

        program = {"enrich": False, "empiric": True, "enrich2": False}
        self.assertEqual(parsers.parse_program(program), "empiric")

        with self.assertRaises(ValueError):
            program = {"enrich": False, "empiric": False, "enrich2": False}
            parsers.parse_program(program)


class TestParseWildTypeSequence(TestCase):
    def test_can_read_from_fasta(self):
        path = os.path.join(TEST_DATA_DIR, "fasta", "lower.fa")
        wtseq = parsers.parse_wt_sequence(path, program="enrich2", non_coding=True)
        expected = (
            "ACAGTTGGATATAGTAGTTTGTACGAGTTGCTTGTGGCTT"
            "CGCCAGCGCATACCAGCATAGTAAAGGCAACGGCCTCTGA"
            "GAGGCTACGATCGTGCCTTGTGGCAAGTCTTCGCTCGCAC"
            "GCCCTTCCTACCGTGCTATGAGAGGAAATCTCGGGCGTA"
        ).upper()
        self.assertEqual(wtseq, expected)

    def test_error_invalid_chars(self):
        with self.assertRaises(exceptions.InvalidWildTypeSequence):
            parsers.parse_wt_sequence("ATXG", program="enrich2", non_coding=True)

    def test_error_not_divisible_by_three_enrich_empiric(self):
        with self.assertRaises(exceptions.SequenceFrameError):
            parsers.parse_wt_sequence("ATGG", program="enrich")
        with self.assertRaises(exceptions.SequenceFrameError):
            parsers.parse_wt_sequence("ATGG", program="empiric")

    def test_error_not_divisible_by_three_enrich2_coding(self):
        with self.assertRaises(exceptions.SequenceFrameError):
            parsers.parse_wt_sequence("ATGG", program="enrich2")

    def test_ok_not_divisible_by_three_enrich2_noncoding(self):
        parsers.parse_wt_sequence("ATGG", program="enrich2", non_coding=True)

    def test_ok_divisible_by_three_enrich2_coding(self):
        parsers.parse_wt_sequence("ATGATC", program="enrich2")

    def test_ok_divisible_by_three_enrich_empiric(self):
        parsers.parse_wt_sequence("ATGATC", program="enrich")
        parsers.parse_wt_sequence("ATGATC", program="empiric")


class TestParseInputType(TestCase):
    @mock.patch("mavedbconvert.parsers.parse_string", return_value="counts")
    def test_calls_parse_string(self, patch):
        parsers.parse_input_type(constants.count_type)
        patch.assert_called()

    def test_error_unrecognised_input_type(self):
        with self.assertRaises(ValueError):
            parsers.parse_input_type("aaa")

    def test_ok_recognised_input_type(self):
        for v in (constants.score_type, constants.count_type):
            parsers.parse_input_type(v)


class TestParseScoreColumn(TestCase):
    @mock.patch("mavedbconvert.parsers.parse_string", return_value="score")
    def test_calls_parse_string(self, patch):
        parsers.parse_score_column("score", constants.score_type, program="enrich")
        patch.assert_called()

    def test_error_enrich_scores_input_and_column_not_defined(self):
        with self.assertRaises(ValueError):
            parsers.parse_score_column(
                value=None, input_type=constants.score_type, program="enrich"
            )

    def test_error_empiric_scores_input_and_column_not_defined(self):
        with self.assertRaises(ValueError):
            parsers.parse_score_column(
                value=None, input_type=constants.score_type, program="empiric"
            )

    def test_ok_enrich_count_input_and_column_not_defined(self):
        parsers.parse_score_column(
            value=None, input_type=constants.count_type, program="enrich"
        )

    def test_ok_empiric_counts_input_and_column_not_defined(self):
        parsers.parse_score_column(
            value=None, input_type=constants.count_type, program="empiric"
        )

    def test_ok_enrich2_and_column_not_defined(self):
        parsers.parse_score_column(
            value=None, input_type=constants.score_type, program="enrich2"
        )


class TestParseOffset(TestCase):
    @mock.patch("mavedbconvert.parsers.parse_numeric", return_value=0)
    def test_calls_parse_numeric(self, patch):
        parsers.parse_offset(0, program="enrich")
        patch.assert_called()

    def test_error_enrich2_is_coding_and_not_mult_of_three(self):
        with self.assertRaises(ValueError):
            parsers.parse_offset(1, "enrich2", non_coding=False)

    def test_ok_enrich2_is_coding_and_mult_of_three(self):
        self.assertEqual(-6, parsers.parse_offset("-6", "enrich2", non_coding=False))

    def test_ok_enrich2_non_coding_and_not_mult_of_three(self):
        self.assertEqual(-7, parsers.parse_offset("-7", "enrich2", non_coding=True))

    def test_error_enrich_empiric_offset_not_mult_of_three(self):
        with self.assertRaises(ValueError):
            parsers.parse_offset(1, "enrich", non_coding=False)

        with self.assertRaises(ValueError):
            parsers.parse_offset(1, "empiric", non_coding=False)

    def test_ok_enrich_empiric_offset_mult_of_three(self):
        self.assertEqual(-6, parsers.parse_offset("-6", "enrich"))
        self.assertEqual(-6, parsers.parse_offset("-6", "empiric"))


class TestParseDocopt(TestCase):
    @staticmethod
    def mock_args(
        program=None,
        src=None,
        dst=None,
        wtseq="AAA",
        offset=0,
        zero_based=False,
        score_column="score",
        hgvs_column=None,
        input_type="scores",
        skip_header="0",
        skip_footer="0",
        sheet_name=None,
        non_coding=False,
    ):

        if program is None:
            program = "enrich2"
        if src is None:
            src = os.path.join(TEST_DATA_DIR, "enrich2", "enrich2.tsv")
        return {
            "enrich": True if program == "enrich" else False,
            "enrich2": True if program == "enrich2" else False,
            "empiric": True if program == "empiric" else False,
            "<src>": os.path.join(TEST_DATA_DIR, program, src),
            "--dst": os.path.join(TEST_DATA_DIR, program, dst) if dst else dst,
            "--score-column": score_column,
            "--hgvs-column": hgvs_column,
            "--skip-header": skip_header,
            "--skip-footer": skip_footer,
            "--sheet-name": sheet_name,
            "--wtseq": wtseq,
            "--offset": offset,
            "--input-type": input_type,
            "--zero-based": zero_based,
            "--non-coding": non_coding,
        }

    def test_returns_correct_program(self):
        for p in ("enrich", "enrich2", "empiric"):
            args = self.mock_args(program=p)
            self.assertEqual(parsers.parse_docopt(args)[0], p)

    def test_is_coding_is_flip_of_non_coding(self):
        args = self.mock_args(non_coding=False)
        _, kwargs = parsers.parse_docopt(args)
        self.assertTrue(kwargs["is_coding"])

        args = self.mock_args(non_coding=True)
        _, kwargs = parsers.parse_docopt(args)
        self.assertFalse(kwargs["is_coding"])

    def test_one_based_is_flip_of_zero_based(self):
        args = self.mock_args(zero_based=False)
        _, kwargs = parsers.parse_docopt(args)
        self.assertTrue(kwargs["one_based"])

        args = self.mock_args(zero_based=True)
        _, kwargs = parsers.parse_docopt(args)
        self.assertFalse(kwargs["one_based"])

    def test_contains_wt_sequence_key(self):
        args = self.mock_args()
        _, kwargs = parsers.parse_docopt(args)
        self.assertIn("wt_sequence", kwargs)

    def test_contains_skip_footer_rows_key(self):
        args = self.mock_args()
        _, kwargs = parsers.parse_docopt(args)
        self.assertIn("skip_footer_rows", kwargs)

    def test_contains_skip_header_rows_key(self):
        args = self.mock_args()
        _, kwargs = parsers.parse_docopt(args)
        self.assertIn("skip_header_rows", kwargs)
