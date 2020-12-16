import os
from unittest.mock import patch
from itertools import product

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from mavedbconvert import validators, enrich2, constants

from tests import ProgramTestCase


class TestEnrich2ParseInput(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.wt = "GCTGAT"
        self.path = os.path.join(self.data_dir, "enrich2", "test_store.h5")
        self.store = pd.HDFStore(self.path, "w")
        self.enrich2 = enrich2.Enrich2(
            self.path, wt_sequence=self.wt, offset=0, one_based=True
        )

        scores, shared, counts, *_ = self.mock_variants_frames()
        self.store["/main/variants/scores/"] = scores
        self.store["/main/variants/scores_shared/"] = shared
        self.store["/main/variants/counts/"] = counts

        scores, shared, counts, *_ = self.mock_synonymous_frames()
        self.store["/main/synonymous/scores/"] = scores
        self.store["/main/synonymous/scores_shared/"] = shared
        self.store["/main/synonymous/counts/"] = counts

        self.files = [
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_synonymous_counts_c1.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_synonymous_counts_c2.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_synonymous_scores_c1.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_synonymous_scores_c2.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_variants_counts_c1.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_variants_counts_c2.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_variants_scores_c1.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_variants_scores_c2.csv",
                )
            ),
        ]

        self.store.close()
        self.store = pd.HDFStore(self.path, mode="r")

    def tearDown(self):
        self.store.close()

    def mock_variants_frames(self, scores_hgvs=None, counts_hgvs=None):
        counts_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["t0", "t1"]],
            names=["condition", "selection", "timepoint"],
        )
        scores_shared_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["SE", "score"]],
            names=["condition", "selection", "value"],
        )
        scores_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["SE", "epsilon", "score"]], names=["condition", "value"]
        )

        if scores_hgvs is None:
            scores_hgvs = [
                "c.2C>T (p.Ala1Val), c.3T>C (p.Ala1=)",
                "c.5A>G (p.Asp2Gly), c.6T>A (p.Asp2Glu)",
            ]
        if counts_hgvs is None:
            counts_hgvs = [
                "c.2C>T (p.Ala1Val), c.3T>C (p.Ala1=)",
                "c.5A>G (p.Asp2Gly), c.6T>A (p.Asp2Glu)",
            ]

        expected = self.parse_rows(scores_hgvs)
        expected_nt = [t[0] for t in expected]
        expected_pro = [t[1] for t in expected]

        scores = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_index)),
            index=scores_hgvs,
            columns=scores_index,
        )
        shared = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_shared_index)),
            index=scores_hgvs,
            columns=scores_shared_index,
        )
        counts = pd.DataFrame(
            np.random.randint(
                low=0, high=100, size=(len(scores_hgvs), len(counts_index))
            ),
            index=counts_hgvs,
            columns=counts_index,
        )
        return scores, shared, counts, expected_nt, expected_pro

    def mock_synonymous_frames(self, scores_hgvs=None, counts_hgvs=None):
        counts_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["t0", "t1"]],
            names=["condition", "selection", "timepoint"],
        )
        scores_shared_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["SE", "score"]],
            names=["condition", "selection", "value"],
        )
        scores_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["SE", "epsilon", "score"]], names=["condition", "value"]
        )

        if scores_hgvs is None:
            scores_hgvs = ["p.Ala1Val, p.Ala1=", "p.Asp2Gly, p.Asp2Glu"]
        if counts_hgvs is None:
            counts_hgvs = ["p.Ala1Val, p.Ala1=", "p.Asp2Gly, p.Asp2Glu"]

        expected = self.parse_rows(scores_hgvs)
        expected_nt = [t[0] for t in expected]
        expected_pro = [t[1] for t in expected]

        scores = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_index)),
            index=scores_hgvs,
            columns=scores_index,
        )
        shared = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_shared_index)),
            index=scores_hgvs,
            columns=scores_shared_index,
        )
        counts = pd.DataFrame(
            np.random.randint(
                low=0, high=100, size=(len(scores_hgvs), len(counts_index))
            ),
            index=counts_hgvs,
            columns=counts_index,
        )
        return scores, shared, counts, expected_nt, expected_pro

    def tearDown(self):
        self.store.close()
        super().tearDown()
        if os.path.isdir(self.enrich2.output_directory):
            os.removedirs(self.enrich2.output_directory)

    def parse_rows(self, variants, element=None):
        return [self.enrich2.parse_row((v, element)) for v in list(variants)]

    @patch.object(pd.DataFrame, "to_csv", return_value=None)
    def test_saves_to_output_directory(self, patch):
        output = os.path.join(self.data_dir, "enrich2", "new")
        p = enrich2.Enrich2(src=self.store, dst=output, wt_sequence=self.wt, offset=0)
        p.parse_input(p.load_input_file())
        for call_args in patch.call_args_list:
            self.assertIn(output, call_args[0][0])

    @patch.object(pd.DataFrame, "to_csv", return_value=None)
    def test_saves_to_file_location_if_no_dst_supplied(self, patch):
        p = enrich2.Enrich2(src=self.store, wt_sequence=self.wt, offset=0)
        p.parse_input(self.enrich2.load_input_file())
        expected_base_path = os.path.normpath(
            os.path.join(self.data_dir, "enrich2", "test_store")
        )
        for call_args in patch.call_args_list:
            self.assertIn(expected_base_path, call_args[0][0])

    @patch("mavedbconvert.enrich2.get_replicate_score_dataframes")
    def test_iterates_over_all_available_tables(self, patch):
        self.enrich2.convert()
        self.assertIn(constants.synonymous_table, patch.call_args_list[0][0])
        self.assertIn(constants.variants_table, patch.call_args_list[1][0])

    @patch(
        "mavedbconvert.enrich2.drop_null",
        side_effect=lambda scores_df, counts_df: (scores_df, counts_df),
    )
    def test_calls_drop_null(self, patch):
        self.enrich2.convert()
        patch.assert_called()

    def test_scores_index_order_retained_in_hgvs_columns(self):
        self.enrich2.convert()

        *_, expected_nt, expected_pro = self.mock_variants_frames()
        nt_pro_tuples = self.parse_rows(
            self.store["/main/variants/scores/"]["c1"].index
        )
        self.assertListEqual(expected_nt, [t[0] for t in nt_pro_tuples])
        self.assertListEqual(expected_pro, [t[1] for t in nt_pro_tuples])

        *_, expected_nt, expected_pro = self.mock_synonymous_frames()
        nt_pro_tuples = self.parse_rows(
            self.store["/main/synonymous/scores/"]["c1"].index
        )
        self.assertListEqual(expected_nt, [t[0] for t in nt_pro_tuples])
        self.assertListEqual(expected_pro, [t[1] for t in nt_pro_tuples])

    def test_counts_index_order_retained_in_hgvs_columns(self):
        self.enrich2.convert()

        *_, expected_nt, expected_pro = self.mock_variants_frames()
        nt_pro_tuples = self.parse_rows(
            self.store["/main/variants/counts/"]["c1"].index
        )
        self.assertListEqual(expected_nt, [t[0] for t in nt_pro_tuples])
        self.assertListEqual(expected_pro, [t[1] for t in nt_pro_tuples])

        *_, expected_nt, expected_pro = self.mock_synonymous_frames()
        nt_pro_tuples = self.parse_rows(
            self.store["/main/synonymous/counts/"]["c1"].index
        )
        self.assertListEqual(expected_nt, [t[0] for t in nt_pro_tuples])
        self.assertListEqual(expected_pro, [t[1] for t in nt_pro_tuples])

    def test_outputs_expected_synonymous_counts_for_each_condition(self):
        self.enrich2.convert()
        *_, _, expected_pro = self.mock_synonymous_frames()

        # C1
        result = pd.read_csv(self.files[0], sep=",")
        expected = pd.DataFrame({constants.pro_variant_col: expected_pro})
        for (rep, tp) in product(["rep1", "rep2"], ["t0", "t1"]):
            expected[rep + "_" + tp] = self.store["/main/synonymous/counts/"]["c1"][
                rep
            ][tp].values.astype(int)
        assert_frame_equal(result, expected)

        # C2
        result = pd.read_csv(self.files[1], sep=",")
        expected = pd.DataFrame({constants.pro_variant_col: expected_pro})
        for (rep, tp) in product(["rep1", "rep2"], ["t0", "t1"]):
            expected[rep + "_" + tp] = self.store["/main/synonymous/counts/"]["c2"][
                rep
            ][tp].values.astype(int)
        assert_frame_equal(result, expected)

    def test_outputs_expected_synonymous_scores_for_each_condition(self):
        self.enrich2.convert()
        *_, _, expected_pro = self.mock_synonymous_frames()
        table_scores = "/main/synonymous/scores/"
        table_shared = "/main/synonymous/scores_shared/"

        # C1
        result = pd.read_csv(self.files[2], sep=",")
        expected = pd.DataFrame(
            {
                constants.pro_variant_col: expected_pro,
                "SE": self.store[table_scores]["c1"]["SE"].values.astype(float),
                "epsilon": self.store[table_scores]["c1"]["epsilon"].values.astype(
                    float
                ),
                "score": self.store[table_scores]["c1"]["score"].values.astype(float),
            },
            columns=[
                constants.pro_variant_col,
                "SE",
                "epsilon",
                "score",
                "SE_rep1",
                "score_rep1",
                "SE_rep2",
                "score_rep2",
            ],
        )
        for (value, rep) in product(["SE", "score"], ["rep1", "rep2"]):
            expected[value + "_" + rep] = self.store[table_shared]["c1"][rep][
                value
            ].values.astype(float)
        assert_frame_equal(result, expected)

        # C2
        result = pd.read_csv(self.files[3], sep=",")
        expected = pd.DataFrame(
            {
                constants.pro_variant_col: expected_pro,
                "SE": self.store[table_scores]["c2"]["SE"].values.astype(float),
                "epsilon": self.store[table_scores]["c2"]["epsilon"].values.astype(
                    float
                ),
                "score": self.store[table_scores]["c2"]["score"].values.astype(float),
            },
            columns=[
                constants.pro_variant_col,
                "SE",
                "epsilon",
                "score",
                "SE_rep1",
                "score_rep1",
                "SE_rep2",
                "score_rep2",
            ],
        )
        for (value, rep) in product(["SE", "score"], ["rep1", "rep2"]):
            expected[value + "_" + rep] = self.store[table_shared]["c2"][rep][
                value
            ].values.astype(float)
        assert_frame_equal(result, expected)

    def test_outputs_expected_variants_counts_for_each_condition(self):
        self.enrich2.convert()
        *_, expected_nt, expected_pro = self.mock_variants_frames()

        # C1
        result = pd.read_csv(self.files[4], sep=",")
        expected = pd.DataFrame(
            {
                constants.nt_variant_col: expected_nt,
                constants.pro_variant_col: expected_pro,
            }
        )
        for (rep, tp) in product(["rep1", "rep2"], ["t0", "t1"]):
            expected[rep + "_" + tp] = self.store["/main/variants/counts/"]["c1"][rep][
                tp
            ].values.astype(int)
        assert_frame_equal(result, expected)

        # C2
        result = pd.read_csv(self.files[5], sep=",")
        expected = pd.DataFrame(
            {
                constants.nt_variant_col: expected_nt,
                constants.pro_variant_col: expected_pro,
            }
        )
        for (rep, tp) in product(["rep1", "rep2"], ["t0", "t1"]):
            expected[rep + "_" + tp] = self.store["/main/variants/counts/"]["c2"][rep][
                tp
            ].values.astype(int)
        assert_frame_equal(result, expected)

    def test_outputs_expected_variants_scores_for_each_condition(self):
        self.enrich2.convert()
        *_, expected_nt, expected_pro = self.mock_variants_frames()
        table_scores = "/main/variants/scores/"
        table_shared = "/main/variants/scores_shared/"

        # C1
        result = pd.read_csv(self.files[6], sep=",")
        expected = pd.DataFrame(
            {
                constants.pro_variant_col: expected_pro,
                constants.nt_variant_col: expected_nt,
                "SE": self.store[table_scores]["c1"]["SE"].values.astype(float),
                "epsilon": self.store[table_scores]["c1"]["epsilon"].values.astype(
                    float
                ),
                "score": self.store[table_scores]["c1"]["score"].values.astype(float),
            },
            columns=[
                constants.nt_variant_col,
                constants.pro_variant_col,
                "SE",
                "epsilon",
                "score",
                "SE_rep1",
                "score_rep1",
                "SE_rep2",
                "score_rep2",
            ],
        )
        for (value, rep) in product(["SE", "score"], ["rep1", "rep2"]):
            expected[value + "_" + rep] = self.store[table_shared]["c1"][rep][
                value
            ].values.astype(float)
        assert_frame_equal(result, expected)

        # C2
        result = pd.read_csv(self.files[7], sep=",")
        expected = pd.DataFrame(
            {
                constants.pro_variant_col: expected_pro,
                constants.nt_variant_col: expected_nt,
                "SE": self.store[table_scores]["c2"]["SE"].values.astype(float),
                "epsilon": self.store[table_scores]["c2"]["epsilon"].values.astype(
                    float
                ),
                "score": self.store[table_scores]["c2"]["score"].values.astype(float),
            },
            columns=[
                constants.nt_variant_col,
                constants.pro_variant_col,
                "SE",
                "epsilon",
                "score",
                "SE_rep1",
                "score_rep1",
                "SE_rep2",
                "score_rep2",
            ],
        )
        for (value, rep) in product(["SE", "score"], ["rep1", "rep2"]):
            expected[value + "_" + rep] = self.store[table_shared]["c2"][rep][
                value
            ].values.astype(float)
        assert_frame_equal(result, expected)

    def test_counts_and_scores_output_define_same_variants_when_input_does_not(self):
        self.store.close()
        self.store = pd.HDFStore(self.path, "w")
        scores, shared, counts, expected_nt, expected_pro = self.mock_variants_frames(
            counts_hgvs=[
                "c.2C>T (p.Ala1Val), c.3T>C (p.Ala1=)",
                # Does not appear in scores
                "c.5A>G (p.Asp2Gly), c.6T>C (p.Asp2=)",
            ]
        )
        self.store["/main/variants/scores/"] = scores
        self.store["/main/variants/scores_shared/"] = shared
        self.store["/main/variants/counts/"] = counts
        self.store.close()
        self.enrich2.convert()

        df_counts = pd.read_csv(self.files[4])  # c1
        df_scores = pd.read_csv(self.files[6])  # c1
        validators.validate_datasets_define_same_variants(df_scores, df_counts)

        df_counts = pd.read_csv(self.files[5])  # c2
        df_scores = pd.read_csv(self.files[7])  # c2
        validators.validate_datasets_define_same_variants(df_scores, df_counts)

    def test_drops_null_rows(self):
        self.store.close()
        self.store = pd.HDFStore(self.path, "w")
        scores, shared, counts, expected_nt, expected_pro = self.mock_variants_frames()

        # Add a null row
        scores = scores.reindex(scores.index.values.tolist() + ["c.1G>G (p.Ala1=)"])
        shared = shared.reindex(shared.index.values.tolist() + ["c.1G>G (p.Ala1=)"])
        counts = counts.reindex(counts.index.values.tolist() + ["c.1G>G (p.Ala1=)"])
        self.store["/main/variants/scores/"] = scores
        self.store["/main/variants/scores_shared/"] = shared
        self.store["/main/variants/counts/"] = counts
        self.store.close()
        self.enrich2.convert()

        df_counts = pd.read_csv(self.files[4])  # c1
        df_scores = pd.read_csv(self.files[6])  # c1
        self.assertNotIn("c.1G>G", df_counts[constants.nt_variant_col])
        self.assertNotIn("c.1G>G", df_scores[constants.nt_variant_col])
        self.assertNotIn("p.Ala1=", df_counts[constants.pro_variant_col])
        self.assertNotIn("p.Ala1=", df_scores[constants.pro_variant_col])

        df_counts = pd.read_csv(self.files[5])  # c1
        df_scores = pd.read_csv(self.files[7])  # c1
        self.assertNotIn("c.1G>G", df_counts[constants.nt_variant_col])
        self.assertNotIn("c.1G>G", df_scores[constants.nt_variant_col])
        self.assertNotIn("p.Ala1=", df_counts[constants.pro_variant_col])
        self.assertNotIn("p.Ala1=", df_scores[constants.pro_variant_col])


class TestEnrich2ParseInputNoVariants(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.wt = "GCTGAT"
        self.path = os.path.join(self.data_dir, "enrich2", "test_store.h5")
        self.store = pd.HDFStore(self.path, "w")
        self.enrich2 = enrich2.Enrich2(
            self.path, wt_sequence=self.wt, offset=0, one_based=True
        )

        scores, shared, counts, *_ = self.mock_synonymous_frames()
        self.store["/main/synonymous/scores/"] = scores
        self.store["/main/synonymous/scores_shared/"] = shared
        self.store["/main/synonymous/counts/"] = counts

        self.files = [
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_synonymous_counts_c1.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_synonymous_counts_c2.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_synonymous_scores_c1.csv",
                )
            ),
            os.path.normpath(
                os.path.join(
                    self.data_dir,
                    "enrich2",
                    "test_store",
                    "mavedb_test_store_synonymous_scores_c2.csv",
                )
            ),
        ]

        self.store.close()
        self.store = pd.HDFStore(self.path, mode="r")

    def mock_synonymous_frames(self, scores_hgvs=None, counts_hgvs=None):
        counts_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["t0", "t1"]],
            names=["condition", "selection", "timepoint"],
        )
        scores_shared_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["rep1", "rep2"], ["SE", "score"]],
            names=["condition", "selection", "value"],
        )
        scores_index = pd.MultiIndex.from_product(
            [["c1", "c2"], ["SE", "epsilon", "score"]], names=["condition", "value"]
        )

        if scores_hgvs is None:
            scores_hgvs = ["p.Ala1Val, p.Ala1=", "p.Asp2Gly, p.Asp2Glu"]
        if counts_hgvs is None:
            counts_hgvs = ["p.Ala1Val, p.Ala1=", "p.Asp2Gly, p.Asp2Glu"]

        expected = self.parse_rows(scores_hgvs)
        expected_nt = [t[0] for t in expected]
        expected_pro = [t[1] for t in expected]

        scores = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_index)),
            index=scores_hgvs,
            columns=scores_index,
        )
        shared = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_shared_index)),
            index=scores_hgvs,
            columns=scores_shared_index,
        )
        counts = pd.DataFrame(
            np.random.randint(
                low=0, high=100, size=(len(scores_hgvs), len(counts_index))
            ),
            index=counts_hgvs,
            columns=counts_index,
        )
        return scores, shared, counts, expected_nt, expected_pro

    def tearDown(self):
        self.store.close()
        super().tearDown()
        if os.path.isdir(self.enrich2.output_directory):
            os.removedirs(self.enrich2.output_directory)

    def parse_rows(self, variants, element=None):
        return [self.enrich2.parse_row((v, element)) for v in list(variants)]

    def test_fails_when_no_variants(self):
        output = os.path.join(self.data_dir, "enrich2", "new")
        p = enrich2.Enrich2(src=self.store, dst=output, wt_sequence=self.wt, offset=0)
        with self.assertRaises(ValueError) as cm:
            p.parse_input(p.load_input_file())
        self.assertEqual(str(cm.exception), "unable to find variants data in HDF5")


class TestEnrich2ParseTsvInput(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.wt = "GCTGAT"
        self.path = os.path.join(self.data_dir, "enrich2", "enrich2.tsv")
        self.enrich2 = enrich2.Enrich2(
            self.path,
            wt_sequence=self.wt,
            offset=0,
            one_based=True,
            hgvs_column="sequence",
            is_coding=False,
        )

    # TODO: this is a totally inadequate set of tests
    # should be incorporated into a refactored more general set of tests
    def test_parses_variants(self):
        result = self.enrich2.convert()
        expected = pd.read_csv(self.path, sep="\t")
        self.assertListEqual(list(result["hgvs_nt"]), list(expected["sequence"]))
