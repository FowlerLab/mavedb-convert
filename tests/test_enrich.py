import os
import unittest

import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal
from fqfa.constants.iupac.protein import AA_CODES

from mavedbconvert import enrich, constants, utilities

from tests import ProgramTestCase


WT = (
    "GACGTTCCACTGCCGGCTGGTTGGGAAATGGCTAAAACTAGTTCTGGTCAGCGTTACTTC"
    "CTGAACCACATCGACCAGACCACCACGTGGCAGGACCCGCGT"
)


class TestEnrichInit(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich", "enrich2.tsv")

    def test_offset_inframe(self):
        enrich.Enrich(src=self.path, wt_sequence="ATC", offset=3)

    def test_error_offset_not_inframe(self):
        with self.assertRaises(ValueError):
            enrich.Enrich(src=self.path, wt_sequence="ATC", offset=1)

    def test_error_noncoding(self):
        with self.assertRaises(ValueError):
            enrich.Enrich(src=self.path, wt_sequence="ATC", is_coding=False)


class TestEnrichParseRow(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich", "enrich.tsv")
        self.enrich = enrich.Enrich(
            src=self.path,
            wt_sequence=WT,
            one_based=False,
            score_column="A",
            input_type=constants.score_type,
        )

    def test_error_pos_len_not_equal_aa_len(self):
        with self.assertRaises(ValueError):
            self.enrich.parse_row("1,2,3,4-L,Y,T")

    # TODO: Uncomment if normalizing variants.
    # def test_converts_qmark_to_xaa_to_single_q(self):
    #     self.assertEqual(self.enrich.parse_row("1-?"), "p.Val2Xaa")
    #     self.assertEqual(
    #         self.enrich.parse_row("0,1-L,?"), "p.[Asp1Leu;Val2Xaa]"
    #     )

    def test_error_no_events(self):
        with self.assertRaises(ValueError):
            self.enrich.parse_row("-")

    def test_error_malformed_seqid(self):
        with self.assertRaises(ValueError):
            self.enrich.parse_row("NA-NA")

    def test_no_changes_uses_equal_event(self):
        self.assertEqual(self.enrich.parse_row("0-D"), "p.Asp1=")

    def test_output_positions_are_1_based(self):
        self.assertEqual(self.enrich.parse_row("0-L"), "p.Asp1Leu")

    def test_parses_correctly_0_based(self):
        protein_seq = utilities.translate_dna(WT, offset=0)
        self.enrich.one_based = False
        aa1_pos = (0 - int(self.enrich.one_based)) + 1
        aa2_pos = (1 - int(self.enrich.one_based)) + 1
        wt_aa1 = AA_CODES[protein_seq[aa1_pos - 1]]
        wt_aa2 = AA_CODES[protein_seq[aa2_pos - 1]]
        expected = "p.[{}]".format(
            ";".join(
                [
                    wt_aa1 + str(aa1_pos) + AA_CODES["L"],
                    wt_aa2 + str(aa2_pos) + AA_CODES["Y"],
                ]
            )
        )
        result = self.enrich.parse_row("0,1-L,Y")
        self.assertEqual(result, expected)

    def test_parses_correctly_1_based(self):
        protein_seq = utilities.translate_dna(WT, offset=0)
        self.enrich.one_based = True
        aa1_pos = (1 - int(self.enrich.one_based)) + 1
        aa2_pos = (2 - int(self.enrich.one_based)) + 1
        wt_aa1 = AA_CODES[protein_seq[aa1_pos - 1]]
        wt_aa2 = AA_CODES[protein_seq[aa2_pos - 1]]
        expected = "p.[{}]".format(
            ";".join(
                [
                    wt_aa1 + str(aa1_pos) + AA_CODES["L"],
                    wt_aa2 + str(aa2_pos) + AA_CODES["Y"],
                ]
            )
        )
        result = self.enrich.parse_row("1,2-L,Y")
        self.assertEqual(result, expected)

    def test_index_error_out_of_bounds(self):
        with self.assertRaises(IndexError):
            self.enrich.parse_row("100-L")

    def test_index_error_negative_index(self):
        with self.assertRaises(IndexError):
            self.enrich.one_based = True
            self.enrich.parse_row("0-L")

    def test_removes_duplicates(self):
        self.assertEqual(self.enrich.parse_row("0,0-L,L"), "p.Asp1Leu")

    def test_uses_three_qmarks(self):
        self.assertEqual(self.enrich.parse_row("0,1-?,?"), "p.[Asp1???;Val2???]")

    def test_applies_offset_divided_by_3(self):
        self.enrich.offset = -3
        self.assertEqual(self.enrich.parse_row("0-D"), "p.Val2Asp")


class TestEnrichParseInput(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich", "enrich.tsv")
        self.enrich = enrich.Enrich(
            src=self.path,
            wt_sequence=WT,
            one_based=False,
            score_column="A",
            input_type=constants.score_type,
        )

    def test_orders_columns(self):
        df = pd.DataFrame({"seqID": ["0,1,2,3-L,Y,T,I"], "A": [1.2], "B": [2.4]})
        result = self.enrich.parse_input(df)
        self.assertEqual(list(result.columns).index(constants.pro_variant_col), 0)
        self.assertEqual(list(result.columns).index(constants.mavedb_score_column), 1)

    def test_removes_hgvs_nt(self):
        df = pd.DataFrame({"seqID": ["0,1,2,3-L,Y,T,I"], "A": [1.2], "B": [2.4]})
        result = self.enrich.parse_input(df)
        self.assertNotIn(constants.nt_variant_col, result.columns)

    def test_removes_null_columns(self):
        df = pd.DataFrame({"seqID": ["0,1,2,3-L,Y,T,I"], "B": [None], "A": [2.4]})
        result = self.enrich.parse_input(df)
        self.assertNotIn("B", result.columns)

    def test_removes_null_rows(self):
        df = pd.DataFrame(
            {"seqID": ["0,1-L,Y", "2,3-T,I"], "A": [None, 1.2], "B": [None, 2.4]}
        )
        result = self.enrich.parse_input(df)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[constants.mavedb_score_column].values[0], 1.2)
        self.assertEqual(result["B"].values[0], 2.4)

    def test_renames_score_column_to_score_and_drops_original(self):
        df = pd.DataFrame({"seqID": ["0,1,2,3-L,Y,T,I"], "A": [1.2], "B": [2.4]})
        result = self.enrich.parse_input(df)
        self.assertListEqual(list(df["A"]), list(result["score"]))
        self.assertIn("B", result.columns)
        self.assertNotIn("A", result.columns)

    def test_keeps_int_type_as_int(self):
        df = pd.DataFrame({"seqID": ["0,1,2,3-L,Y,T,I"], "A": [1]})
        result = self.enrich.parse_input(df)
        self.assertTrue(
            np.issubdtype(
                result[constants.mavedb_score_column].values[0], np.signedinteger
            )
        )

    def test_removes_non_numeric(self):
        df = pd.DataFrame({"seqID": ["0,1,2,3-L,Y,T,I"], "A": [1.2], "B": ["a"]})
        result = self.enrich.parse_input(df)
        self.assertNotIn("B", result)


class TestEnrichLoadInput(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich", "enrich.tsv")
        self.path_1based = os.path.join(self.data_dir, "enrich", "enrich_1based.tsv")
        self.path_csv = os.path.join(self.data_dir, "enrich", "enrich1.csv")
        self.expected = os.path.join(self.data_dir, "enrich", "enrich_expected.csv")
        self.expected_offset = os.path.join(
            self.data_dir, "enrich", "enrich_expected_offset.csv"
        )
        self.excel_path = os.path.join(self.data_dir, "enrich", "enrich.xlsx")
        self.no_seq_id = os.path.join(self.data_dir, "enrich", "enrich_no_seqid.tsv")
        self.tmp_path = os.path.join(self.data_dir, "enrich", "tmp.xlsx")

    def test_error_seq_id_not_in_columns(self):
        p = enrich.Enrich(
            src=self.no_seq_id,
            wt_sequence=WT,
            score_column="log2_ratio",
            input_type=constants.score_type,
        )
        with self.assertRaises(ValueError):
            p.load_input_file()

    def test_loads_first_sheet_by_default(self):
        data = [
            {"seqID": ["1,2,3-G,L-Y"], "score": [1.2]},
            {"seqID": ["1,2,5-G,L-Y"], "score": [1.4]},
        ]
        self.mock_multi_sheet_excel_file(self.tmp_path, data)
        p = enrich.Enrich(
            src=self.tmp_path,
            wt_sequence=WT,
            score_column="score",
            input_type=constants.score_type,
        )
        df = p.load_input_file()
        expected = pd.DataFrame(data[0])
        assert_frame_equal(df, expected)

    def test_loads_xlxs(self):
        p = enrich.Enrich(
            src=self.excel_path,
            wt_sequence=WT,
            score_column="log2_ratio",
            input_type=constants.score_type,
        )
        result = p.load_input_file()
        expected = pd.read_excel(self.excel_path, na_values=constants.extra_na)
        assert_frame_equal(result, expected)

    def test_loads_table(self):
        p = enrich.Enrich(
            src=self.path,
            wt_sequence=WT,
            score_column="log2_ratio",
            input_type=constants.score_type,
        )
        result = p.load_input_file()
        expected = pd.read_csv(self.path, delimiter="\t", na_values=constants.extra_na)
        assert_frame_equal(result, expected)

    def test_loads_csv(self):
        expected = pd.read_csv(self.path, delimiter="\t", na_values=constants.extra_na)
        expected.to_csv(self.path_csv, index=False)
        p = enrich.Enrich(
            src=self.path_csv,
            wt_sequence=WT,
            score_column="log2_ratio",
            input_type=constants.score_type,
        )
        result = p.load_input_file()
        assert_frame_equal(result, expected)

    def test_table_and_excel_load_same_dataframe(self):
        p1 = enrich.Enrich(
            src=self.path,
            wt_sequence=WT,
            score_column="log2_ratio",
            input_type=constants.score_type,
        )
        p2 = enrich.Enrich(
            src=self.excel_path,
            wt_sequence=WT,
            score_column="log2_ratio",
            input_type=constants.score_type,
        )
        assert_frame_equal(p1.load_input_file(), p2.load_input_file())


class TestEnrichIntegration(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(self.data_dir, "enrich", "enrich.tsv")
        self.path_1based = os.path.join(self.data_dir, "enrich", "enrich_1based.tsv")
        self.excel_path = os.path.join(self.data_dir, "enrich", "enrich.xlsx")
        self.no_seq_id = os.path.join(self.data_dir, "enrich", "enrich_no_seqid.tsv")

        self.expected = os.path.join(self.data_dir, "enrich", "enrich_expected.csv")
        self.expected_offset = os.path.join(
            self.data_dir, "enrich", "enrich_expected_offset.csv"
        )

    def test_saves_to_input_dst_by_default(self):
        p = enrich.Enrich(
            src=self.path,
            wt_sequence=WT,
            one_based=False,
            score_column="log2_ratio",
            input_type=constants.score_type,
        )
        p.convert()
        self.assertTrue(
            os.path.isfile(os.path.join(self.data_dir, "enrich", "mavedb_enrich.csv"))
        )

    def test_output_with_offset(self):
        p = enrich.Enrich(
            src=self.path,
            wt_sequence=WT,
            offset=-3,
            one_based=False,
            score_column="log2_ratio",
            input_type=constants.score_type,
        )
        p.convert()
        result = pd.read_csv(os.path.join(self.data_dir, "enrich", "mavedb_enrich.csv"))
        expected = pd.read_csv(self.expected_offset)
        assert_frame_equal(expected, result)

    def test_output_from_one_based_input(self):
        p = enrich.Enrich(
            src=self.path_1based,
            wt_sequence=WT,
            one_based=True,
            score_column="log2_ratio",
            input_type=constants.score_type,
        )
        p.convert()
        result = pd.read_csv(
            os.path.join(self.data_dir, "enrich", "mavedb_enrich_1based.csv")
        )
        expected = pd.read_csv(self.expected)
        assert_frame_equal(expected, result)


if __name__ == "__main__":
    unittest.main()
