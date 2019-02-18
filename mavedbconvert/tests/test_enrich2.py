import os
import mock
from unittest import TestCase
from itertools import product

import hgvsp
import hgvsp.constants

import numpy as np
import pandas as pd
from pandas.testing import assert_index_equal, assert_frame_equal

from .. import validators, enrich2, constants, exceptions

from . import ProgramTestCase


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(BASE_DIR + '/data/')


# Utility tests
# --------------------------------------------------------------------------- #
class TestGetCountDataFrames(TestCase):
    """
    Test method get_count_dataframes checking if conditions are correctly
    parsed.
    """
    def setUp(self):
        self.path = os.path.join(DATA_DIR, 'test_store.h5')
        self.store = pd.HDFStore(self.path, 'w')
        index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['rep1', 'rep2'], ['t0', 't1']],
            names=['condition', 'selection', 'timepoint']
        )
        scores_hgvs = ['c.1A>G', 'c.3A>G']
        counts_hgvs = ['c.1A>G', 'c.2A>G']
        self.store['/main/variants/scores/'] = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(index)), index=scores_hgvs,
            columns=index
        )
        self.store['/main/variants/counts/'] = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(index)), index=counts_hgvs,
            columns=index
        )
        
    def tearDown(self):
        self.store.close()
        if os.path.isfile(self.path):
            os.unlink(self.path)
        
    def test_column_names_combine_selection_and_timepoint(self):
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd='c1')
        self.assertListEqual(
            list(cnd_df.columns),
            ['rep1_t0', 'rep1_t1', 'rep2_t0', 'rep2_t1'],
        )

    def test_index_of_dfs_match_index_of_scores(self):
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd='c1')
        assert_index_equal(
            self.store['/main/variants/scores/'].index, cnd_df.index)

    def test_row_filled_with_nans_filtered_index_not_in_counts(self):
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd='c1')
        self.assertTrue(np.all(cnd_df.loc['c.3A>G', :].isnull()))

    def test_returns_empty_when_missing_scores_key(self):
        self.store.remove('/main/variants/scores/')
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd='c1')
        self.assertIsNone(cnd_df)

    def test_returns_empty_when_missing_counts_key(self):
        self.store.remove('/main/variants/counts/')
        cnd_df = enrich2.get_count_dataframe_by_condition(self.store, cnd='c1')
        self.assertIsNone(cnd_df)
        
        
class TestFlattenColumnNames(TestCase):
    def setUp(self):
        index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['rep1', 'rep2'], ['t0', 't1']],
            names=['condition', 'selection', 'timepoint']
        )
        scores_hgvs = ['c.1A>G', 'c.3A>G']
        self.df = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(index)), index=scores_hgvs,
            columns=index
        )
    
    def test_column_names_combine_columns_using_ordering(self):
        cnames = enrich2.flatten_column_names(
            self.df.loc[:, pd.IndexSlice['c1', :, :]].columns, ordering=(2, 1))
        self.assertListEqual(
            cnames, ['t0_rep1', 't1_rep1', 't0_rep2', 't1_rep2'])


class TestReplicateScoreDataFrames(TestCase):
    """
    Test method get_replicate_score_dataframes checking if conditions are
    correctly parsed.
    """
    def setUp(self):
        self.path = os.path.join(DATA_DIR, 'test_store.h5')
        self.store = pd.HDFStore(self.path, 'w')
        
        shared_index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['rep1', 'rep2'], ['SE', 'score']],
            names=['condition', 'selection', 'value']
        )
        index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['SE', 'epsilon', 'score']],
            names=['condition', 'value']
        )
        
        hgvs = ['c.1A>G', 'c.2A>G']
        self.store['/main/variants/scores/'] = pd.DataFrame(
            np.random.randn(len(hgvs), len(index)),
            index=hgvs,
            columns=index
        )
        self.store['/main/variants/scores_shared/'] = pd.DataFrame(
            np.random.randn(len(hgvs), len(shared_index)),
            index=hgvs,
            columns=shared_index
        )

    def tearDown(self):
        self.store.close()
        if os.path.isfile(self.path):
            os.unlink(self.path)
            
    def test_conditions_are_dictionary_keys(self):
        cnd_dfs = enrich2.get_replicate_score_dataframes(self.store)
        self.assertIn('c1', cnd_dfs)
        self.assertIn('c2', cnd_dfs)

    def test_returns_empty_when_missing_scores_key(self):
        self.store.remove('/main/variants/scores')
        cnd_dfs = enrich2.get_replicate_score_dataframes(self.store)
        self.assertDictEqual(cnd_dfs, {})

    def test_returns_empty_when_missing_shared_scores_key(self):
        self.store.remove('/main/variants/scores_shared')
        cnd_dfs = enrich2.get_replicate_score_dataframes(self.store)
        self.assertDictEqual(cnd_dfs, {})

    def test_adds_rep_id_to_score_and_SE(self):
        cnd_dfs = enrich2.get_replicate_score_dataframes(self.store)
        for cnd, df in cnd_dfs.items():
            for c_name in df.columns:
                if c_name not in ('score', 'SE', 'epsilon'):
                    self.assertIn('rep', c_name.lower())

    def test_assertion_error_scores_shared_scores_different_index(self):
        hgvs = ['c.1A>G', 'c.3A>G']
        shared_index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['rep1', 'rep2'], ['SE', 'score']],
            names=['condition', 'selection', 'value']
        )
        self.store['/main/variants/scores_shared/'] = pd.DataFrame(
            np.random.randn(len(hgvs), len(shared_index)),
            index=hgvs,
            columns=shared_index
        )
        with self.assertRaises(AssertionError):
            enrich2.get_replicate_score_dataframes(self.store)


class TestDropNull(TestCase):
    def test_calls_drop_na_rows_from_scores_inplace(self):
        df = pd.DataFrame({'A': [None, 1]})
        enrich2.drop_null(df)
        assert_index_equal(df.index, pd.Index([1]))
        
    def test_calls_drop_na_cols_from_scores_inplace(self):
        df = pd.DataFrame({'A': [1, 2], 'B': [None, None]})
        enrich2.drop_null(df)
        self.assertNotIn('B', df)
        
    def test_assertion_error_if_counts_scores_indices_do_not_match(self):
        df1 = pd.DataFrame({'A': [1, 2]}, index=['a', 'b'])
        df2 = pd.DataFrame({'A': [1, 2]}, index=['c', 'b'])
        with self.assertRaises(AssertionError):
            enrich2.drop_null(df1, df2)
            
    def test_assertion_error_if_dfs_define_different_variants(self):
        df1 = pd.DataFrame(
            {constants.nt_variant_col: ['c.1A>G'],
             constants.pro_variant_col: ['p.G4L']},
            index=['c.1A>G'])
        df2 = pd.DataFrame(
            {constants.nt_variant_col: ['c.2A>G'],
             constants.pro_variant_col: ['p.G4L']},
            index=['c.1A>G'])
        with self.assertRaises(AssertionError):
            enrich2.drop_null(df1, df2)
        
    def test_na_rows_dropped_from_scores_counts_after_join(self):
        df1 = pd.DataFrame({
            constants.nt_variant_col: ['c.1A>G', 'c.2A>G'],
            constants.pro_variant_col: ['p.G4L', 'p.G5L'],
            'score': [1, None]}, index=['c.1A>G', 'c.2A>G']
        )
        df2 = pd.DataFrame({
            constants.nt_variant_col: ['c.1A>G', 'c.2A>G'],
            constants.pro_variant_col: ['p.G4L', 'p.G5L'],
            'count': [10, None]}, index=['c.1A>G', 'c.2A>G']
        )
        scores, counts = enrich2.drop_null(df1, df2)
        self.assertNotIn('c.2A>G', scores.index.values)
        self.assertNotIn('c.2A>G', counts.index.values)
        
    def test_na_cols_dropped_from_scores_counts_after_join(self):
        df1 = pd.DataFrame({
            constants.nt_variant_col: ['c.1A>G', ],
            constants.pro_variant_col: [None, ],
            'score': [1, ]}, index=['c.1A>G', ]
        )
        df2 = pd.DataFrame({
            constants.nt_variant_col: ['c.1A>G', ],
            constants.pro_variant_col: [None, ],
            'count': [10, ]}, index=['c.1A>G', ]
        )
        scores, counts = enrich2.drop_null(df1, df2)
        self.assertNotIn(constants.pro_variant_col, scores.columns)
        self.assertNotIn(constants.pro_variant_col, counts.columns)
        
    def test_scores_and_counts_columns_separated_after_join(self):
        df1 = pd.DataFrame({
            constants.nt_variant_col: ['c.1A>G', 'c.2A>G'],
            constants.pro_variant_col: ['p.G4L', 'p.G5L'],
            'score': [1, None]}, index=['c.1A>G', 'c.2A>G']
        )
        df2 = pd.DataFrame({
            constants.nt_variant_col: ['c.1A>G', 'c.2A>G'],
            constants.pro_variant_col: ['p.G4L', 'p.G5L'],
            'count': [10, None]}, index=['c.1A>G', 'c.2A>G']
        )
        scores, counts = enrich2.drop_null(df1, df2)
        self.assertListEqual(
            list(scores.columns),
            [constants.nt_variant_col, constants.pro_variant_col, 'score'])
        self.assertListEqual(
            list(counts.columns),
            [constants.nt_variant_col, constants.pro_variant_col, 'count'])
        

# HD5/Row parsing tests
# --------------------------------------------------------------------------- #
class TestEnrich2ConvertH5Filepath(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(DATA_DIR, 'enrich2.h5')
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence='AAA')
        self.bin.append(self.path.replace('.h5', ''))

    def test_replaces_underscore_with_spaces(self):
        res = self.enrich2.convert_h5_filepath(
            'base', 'syn vars', 'scores', 'c1')
        self.assertIn('syn_vars', res)

    def test_concats_basename_elem_type_then_cnd_and_csv_ext(self):
        res = self.enrich2.convert_h5_filepath(
            'base', 'syn vars', 'scores', 'c1')
        expected = 'mavedb_base_syn_vars_scores_c1.csv'
        self.assertEqual(
            res,
            os.path.join(self.enrich2.output_directory, expected)
        )


class TestEnrich2ConvertH5Df(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(DATA_DIR, 'enrich2.h5')
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence='AAA')

    def test_drops_non_numeric_columns(self):
        df = pd.DataFrame(
            data={'score': [1, ], 'B': ['a', ], },
            index=['c.1A>G (p.Lys1Val)']
        )
        result = self.enrich2.convert_h5_df(
            df=df,
            element=constants.variants_table,
            df_type=constants.score_type
        )
        self.assertNotIn('B', result)

    def test_type_casts_numeric_to_int_and_float(self):
        df = pd.DataFrame(
            data={'score': [1, ], 'B': [1.2, ]},
            index=['c.1A>G (p.Lys1Val)'],
        )
        result = self.enrich2.convert_h5_df(
            df=df,
            element=constants.variants_table,
            df_type=constants.score_type
        )
        self.assertTrue(np.issubdtype(
            result['score'].values[0], np.signedinteger))
        self.assertTrue(np.issubdtype(
            result['B'].values[0], np.floating))

    def test_sets_index_as_input_index(self):
        df = pd.DataFrame({
            'score': [1, ], 'B': ['a', ]},
            index=['c.1A>G (p.Lys1Val)']
        )
        result = self.enrich2.convert_h5_df(
            df=df,
            element=constants.variants_table,
            df_type=constants.score_type
        )
        assert_index_equal(result.index, df.index)

    
class TestEnrich2ParseInput(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.wt = 'GCTGAT'
        self.path = os.path.join(DATA_DIR, 'test_store.h5')
        self.store = pd.HDFStore(self.path, 'w')
        self.enrich2 = enrich2.Enrich2(
            self.path, wt_sequence=self.wt, offset=0, one_based=True)

        scores, shared, counts, *_ = self.mock_variants_frames()
        self.store['/main/variants/scores/'] = scores
        self.store['/main/variants/scores_shared/'] = shared
        self.store['/main/variants/counts/'] = counts

        scores, shared, counts, *_ = self.mock_synonymous_frames()
        self.store['/main/synonymous/scores/'] = scores
        self.store['/main/synonymous/scores_shared/'] = shared
        self.store['/main/synonymous/counts/'] = counts

        self.files = [
            os.path.normpath(os.path.join(
                DATA_DIR, 'test_store',
                'mavedb_test_store_synonymous_counts_c1.csv')),
            os.path.normpath(os.path.join(
                DATA_DIR, 'test_store',
                'mavedb_test_store_synonymous_counts_c2.csv')),
            os.path.normpath(os.path.join(
                DATA_DIR, 'test_store',
                'mavedb_test_store_synonymous_scores_c1.csv')),
            os.path.normpath(os.path.join(
                DATA_DIR, 'test_store',
                'mavedb_test_store_synonymous_scores_c2.csv')),

            os.path.normpath(os.path.join(
                DATA_DIR, 'test_store',
                'mavedb_test_store_variants_counts_c1.csv')),
            os.path.normpath(os.path.join(
                DATA_DIR, 'test_store',
                'mavedb_test_store_variants_counts_c2.csv')),
            os.path.normpath(os.path.join(
                DATA_DIR, 'test_store',
                'mavedb_test_store_variants_scores_c1.csv')),
            os.path.normpath(os.path.join(
                DATA_DIR, 'test_store',
                'mavedb_test_store_variants_scores_c2.csv')),
        ]

        self.bin.extend(self.files)
        self.bin.append(self.path)
        self.store.close()
        self.store = pd.HDFStore(self.path, mode='r')

    def mock_variants_frames(self, scores_hgvs=None, counts_hgvs=None):
        counts_index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['rep1', 'rep2'], ['t0', 't1']],
            names=['condition', 'selection', 'timepoint']
        )
        scores_shared_index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['rep1', 'rep2'], ['SE', 'score']],
            names=['condition', 'selection', 'value']
        )
        scores_index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['SE', 'epsilon', 'score']],
            names=['condition', 'value']
        )

        if scores_hgvs is None:
            scores_hgvs = [
                'c.2C>T (p.Ala1Val), c.3T>C (p.Ala1=)',
                'c.5A>G (p.Asp2Gly), c.6T>A (p.Asp2Glu)'
            ]
        if counts_hgvs is None:
            counts_hgvs = [
                'c.2C>T (p.Ala1Val), c.3T>C (p.Ala1=)',
                'c.5A>G (p.Asp2Gly), c.6T>A (p.Asp2Glu)'
            ]

        expected = self.parse_rows(scores_hgvs)
        expected_nt = [t[0] for t in expected]
        expected_pro = [t[1] for t in expected]
        
        scores = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_index)),
            index=scores_hgvs,
            columns=scores_index
        )
        shared = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_shared_index)),
            index=scores_hgvs,
            columns=scores_shared_index
        )
        counts = pd.DataFrame(
            np.random.randint(
                low=0, high=100, size=(len(scores_hgvs), len(counts_index))),
            index=counts_hgvs,
            columns=counts_index
        )
        return scores, shared, counts, expected_nt, expected_pro

    def mock_synonymous_frames(self, scores_hgvs=None, counts_hgvs=None):
        counts_index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['rep1', 'rep2'], ['t0', 't1']],
            names=['condition', 'selection', 'timepoint']
        )
        scores_shared_index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['rep1', 'rep2'], ['SE', 'score']],
            names=['condition', 'selection', 'value']
        )
        scores_index = pd.MultiIndex.from_product(
            [['c1', 'c2'], ['SE', 'epsilon', 'score']],
            names=['condition', 'value']
        )
        
        if scores_hgvs is None:
            scores_hgvs = [
                'p.Ala1Val, p.Ala1=',
                'p.Asp2Gly, p.Asp2Glu',
            ]
        if counts_hgvs is None:
            counts_hgvs = [
                'p.Ala1Val, p.Ala1=',
                'p.Asp2Gly, p.Asp2Glu',
            ]

        expected = self.parse_rows(scores_hgvs)
        expected_nt = [t[0] for t in expected]
        expected_pro = [t[1] for t in expected]
    
        scores = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_index)),
            index=scores_hgvs,
            columns=scores_index
        )
        shared = pd.DataFrame(
            np.random.randn(len(scores_hgvs), len(scores_shared_index)),
            index=scores_hgvs,
            columns=scores_shared_index
        )
        counts = pd.DataFrame(
            np.random.randint(
                low=0, high=100, size=(len(scores_hgvs), len(counts_index))),
            index=counts_hgvs,
            columns=counts_index
        )
        return scores, shared, counts, expected_nt, expected_pro

    def tearDown(self):
        self.store.close()
        super().tearDown()
        if os.path.isdir(self.enrich2.output_directory):
            os.removedirs(self.enrich2.output_directory)

    def parse_rows(self, variants, element=None):
        return [
            self.enrich2.parse_row((v, element))
            for v in list(variants)
        ]

    @mock.patch.object(pd.DataFrame, 'to_csv', return_value=None)
    def test_saves_to_output_directory(self, patch):
        output = os.path.join(DATA_DIR, 'new')
        p = enrich2.Enrich2(
            src=self.store, dst=output, wt_sequence=self.wt, offset=0)
        p.parse_input(p.load_input_file())
        for call_args in patch.call_args_list:
            self.assertIn(output, call_args[0][0])
        self.bin.append(output)

    @mock.patch.object(pd.DataFrame, 'to_csv', return_value=None)
    def test_saves_to_file_location_if_no_dst_supplied(self, patch):
        p = enrich2.Enrich2(src=self.store, wt_sequence=self.wt, offset=0)
        p.parse_input(self.enrich2.load_input_file())
        expected_base_path = os.path.normpath(
            os.path.join(DATA_DIR, 'test_store'))
        for call_args in patch.call_args_list:
            self.assertIn(expected_base_path, call_args[0][0])
            
    @mock.patch('mavedbconvert.enrich2.get_replicate_score_dataframes')
    def test_iterates_over_all_available_tables(self, patch):
        self.enrich2.parse_input(self.enrich2.load_input_file())
        self.assertIn(constants.synonymous_table, patch.call_args_list[0][0])
        self.assertIn(constants.variants_table, patch.call_args_list[1][0])

    @mock.patch('mavedbconvert.enrich2.drop_null',
                side_effect=lambda scores_df, counts_df: (scores_df, counts_df))
    def test_calls_drop_null(self, patch):
        self.enrich2.parse_input(self.enrich2.load_input_file())
        patch.assert_called()

    def test_scores_index_order_retained_in_hgvs_columns(self):
        self.enrich2.parse_input(self.enrich2.load_input_file())
        
        *_, expected_nt, expected_pro = self.mock_variants_frames()
        nt_pro_tuples = self.parse_rows(
            self.store['/main/variants/scores/']['c1'].index)
        self.assertListEqual(expected_nt, [t[0] for t in nt_pro_tuples])
        self.assertListEqual(expected_pro, [t[1] for t in nt_pro_tuples])
        
        *_, expected_nt, expected_pro = self.mock_synonymous_frames()
        nt_pro_tuples = self.parse_rows(
            self.store['/main/synonymous/scores/']['c1'].index)
        self.assertListEqual(expected_nt, [t[0] for t in nt_pro_tuples])
        self.assertListEqual(expected_pro, [t[1] for t in nt_pro_tuples])

    def test_counts_index_order_retained_in_hgvs_columns(self):
        self.enrich2.parse_input(self.enrich2.load_input_file())
    
        *_, expected_nt, expected_pro = self.mock_variants_frames()
        nt_pro_tuples = self.parse_rows(
            self.store['/main/variants/counts/']['c1'].index)
        self.assertListEqual(expected_nt, [t[0] for t in nt_pro_tuples])
        self.assertListEqual(expected_pro, [t[1] for t in nt_pro_tuples])
    
        *_, expected_nt, expected_pro = self.mock_synonymous_frames()
        nt_pro_tuples = self.parse_rows(
            self.store['/main/synonymous/counts/']['c1'].index)
        self.assertListEqual(expected_nt, [t[0] for t in nt_pro_tuples])
        self.assertListEqual(expected_pro, [t[1] for t in nt_pro_tuples])
        
    def test_outputs_expected_synonymous_counts_for_each_condition(self):
        self.enrich2.parse_input(self.enrich2.load_input_file())
        *_, _, expected_pro = self.mock_synonymous_frames()
        
        # C1
        result = pd.read_csv(self.files[0], sep=',')
        expected = pd.DataFrame({constants.pro_variant_col: expected_pro, })
        for (rep, tp) in product(['rep1', 'rep2'], ['t0', 't1']):
            expected[rep + '_' + tp] = \
                self.store['/main/synonymous/counts/']['c1'][rep][tp].\
                    values.astype(int)
        assert_frame_equal(result, expected)
        
        # C2
        result = pd.read_csv(self.files[1], sep=',')
        expected = pd.DataFrame({constants.pro_variant_col: expected_pro, })
        for (rep, tp) in product(['rep1', 'rep2'], ['t0', 't1']):
            expected[rep + '_' + tp] = \
                self.store['/main/synonymous/counts/']['c2'][rep][tp].\
                    values.astype(int)
        assert_frame_equal(result, expected)

    def test_outputs_expected_synonymous_scores_for_each_condition(self):
        self.enrich2.parse_input(self.enrich2.load_input_file())
        *_, _, expected_pro = self.mock_synonymous_frames()
        table_scores = '/main/synonymous/scores/'
        table_shared = '/main/synonymous/scores_shared/'
        
        # C1
        result = pd.read_csv(self.files[2], sep=',')
        expected = pd.DataFrame({
            constants.pro_variant_col: expected_pro,
            'SE': self.store[table_scores]['c1']['SE'].values.astype(float),
            'epsilon': self.store[table_scores]['c1']['epsilon'].values.astype(float),
            'score': self.store[table_scores]['c1']['score'].values.astype(float),
        }, columns=[
            constants.pro_variant_col, 'SE', 'epsilon', 'score',
            'SE_rep1', 'score_rep1', 'SE_rep2', 'score_rep2'
        ])
        for (value, rep) in product(['SE', 'score'], ['rep1', 'rep2']):
            expected[value + '_' + rep] = \
                self.store[table_shared]['c1'][rep][value].values.astype(float)
        assert_frame_equal(result, expected)

        # C2
        result = pd.read_csv(self.files[3], sep=',')
        expected = pd.DataFrame({
            constants.pro_variant_col: expected_pro,
            'SE': self.store[table_scores]['c2']['SE'].values.astype(float),
            'epsilon': self.store[table_scores]['c2']['epsilon'].values.astype(float),
            'score': self.store[table_scores]['c2']['score'].values.astype(float),
        }, columns=[
            constants.pro_variant_col, 'SE', 'epsilon', 'score',
            'SE_rep1', 'score_rep1', 'SE_rep2', 'score_rep2'
        ])
        for (value, rep) in product(['SE', 'score'], ['rep1', 'rep2']):
            expected[value + '_' + rep] = \
                self.store[table_shared]['c2'][rep][value].values.astype(float)
        assert_frame_equal(result, expected)

    def test_outputs_expected_variants_counts_for_each_condition(self):
        self.enrich2.parse_input(self.enrich2.load_input_file())
        *_, expected_nt, expected_pro = self.mock_variants_frames()
        
        # C1
        result = pd.read_csv(self.files[4], sep=',')
        expected = pd.DataFrame(
            {constants.nt_variant_col: expected_nt,
             constants.pro_variant_col: expected_pro, }
        )
        for (rep, tp) in product(['rep1', 'rep2'], ['t0', 't1']):
            expected[rep + '_' + tp] = \
                self.store['/main/variants/counts/']['c1'][rep][tp].\
                    values.astype(int)
        assert_frame_equal(result, expected)

        # C2
        result = pd.read_csv(self.files[5], sep=',')
        expected = pd.DataFrame(
            {constants.nt_variant_col: expected_nt,
             constants.pro_variant_col: expected_pro, }
        )
        for (rep, tp) in product(['rep1', 'rep2'], ['t0', 't1']):
            expected[rep + '_' + tp] = \
                self.store['/main/variants/counts/']['c2'][rep][tp].\
                    values.astype(int)
        assert_frame_equal(result, expected)

    def test_outputs_expected_variants_scores_for_each_condition(self):
        self.enrich2.parse_input(self.enrich2.load_input_file())
        *_, expected_nt, expected_pro = self.mock_variants_frames()
        table_scores = '/main/variants/scores/'
        table_shared = '/main/variants/scores_shared/'
    
        # C1
        result = pd.read_csv(self.files[6], sep=',')
        expected = pd.DataFrame({
            constants.pro_variant_col: expected_pro,
            constants.nt_variant_col: expected_nt,
            'SE': self.store[table_scores]['c1']['SE'].values.astype(float),
            'epsilon': self.store[table_scores]['c1']['epsilon'].values.astype(
                float),
            'score': self.store[table_scores]['c1']['score'].values.astype(
                float),
        }, columns=[
            constants.nt_variant_col,
            constants.pro_variant_col, 'SE', 'epsilon', 'score',
            'SE_rep1', 'score_rep1', 'SE_rep2', 'score_rep2'
        ])
        for (value, rep) in product(['SE', 'score'], ['rep1', 'rep2']):
            expected[value + '_' + rep] = \
                self.store[table_shared]['c1'][rep][value].values.astype(float)
        assert_frame_equal(result, expected)
    
        # C2
        result = pd.read_csv(self.files[7], sep=',')
        expected = pd.DataFrame({
            constants.pro_variant_col: expected_pro,
            constants.nt_variant_col: expected_nt,
            'SE': self.store[table_scores]['c2']['SE'].values.astype(float),
            'epsilon': self.store[table_scores]['c2']['epsilon'].values.astype(
                float),
            'score': self.store[table_scores]['c2']['score'].values.astype(
                float),
        }, columns=[
            constants.nt_variant_col,
            constants.pro_variant_col, 'SE', 'epsilon', 'score',
            'SE_rep1', 'score_rep1', 'SE_rep2', 'score_rep2'
        ])
        for (value, rep) in product(['SE', 'score'], ['rep1', 'rep2']):
            expected[value + '_' + rep] = \
                self.store[table_shared]['c2'][rep][value].values.astype(float)
        assert_frame_equal(result, expected)
        
    def test_counts_and_scores_output_define_same_variants_when_input_does_not(self):
        self.store.close()
        self.store = pd.HDFStore(self.path, 'w')
        scores, shared, counts, expected_nt, expected_pro = \
            self.mock_variants_frames(
                counts_hgvs=[
                    'c.2C>T (p.Ala1Val), c.3T>C (p.Ala1=)',
                    # Does not appear in scores
                    'c.5A>G (p.Asp2Gly), c.6T>C (p.Asp2=)',
                ]
            )
        self.store['/main/variants/scores/'] = scores
        self.store['/main/variants/scores_shared/'] = shared
        self.store['/main/variants/counts/'] = counts
        self.store.close()
        self.enrich2.parse_input(self.enrich2.load_input_file())

        df_counts = pd.read_csv(self.files[4])  # c1
        df_scores = pd.read_csv(self.files[6])  # c1
        validators.validate_datasets_define_same_variants(df_scores, df_counts)

        df_counts = pd.read_csv(self.files[5])  # c2
        df_scores = pd.read_csv(self.files[7])  # c2
        validators.validate_datasets_define_same_variants(df_scores, df_counts)
        
    def test_drops_null_rows(self):
        self.store.close()
        self.store = pd.HDFStore(self.path, 'w')
        scores, shared, counts, expected_nt, expected_pro = \
            self.mock_variants_frames()
        
        # Add a null row
        scores = scores.reindex(scores.index.values.tolist() +
                                ['c.1G>G (p.Ala1=)'])
        shared = shared.reindex(shared.index.values.tolist() +
                                ['c.1G>G (p.Ala1=)'])
        counts = counts.reindex(counts.index.values.tolist() +
                                ['c.1G>G (p.Ala1=)'])
        self.store['/main/variants/scores/'] = scores
        self.store['/main/variants/scores_shared/'] = shared
        self.store['/main/variants/counts/'] = counts
        self.store.close()
        self.enrich2.parse_input(self.enrich2.load_input_file())
    
        df_counts = pd.read_csv(self.files[4])  # c1
        df_scores = pd.read_csv(self.files[6])  # c1
        self.assertNotIn('c.1G>G', df_counts[constants.nt_variant_col])
        self.assertNotIn('c.1G>G', df_scores[constants.nt_variant_col])
        self.assertNotIn('p.Ala1=', df_counts[constants.pro_variant_col])
        self.assertNotIn('p.Ala1=', df_scores[constants.pro_variant_col])
    
        df_counts = pd.read_csv(self.files[5])  # c1
        df_scores = pd.read_csv(self.files[7])  # c1
        self.assertNotIn('c.1G>G', df_counts[constants.nt_variant_col])
        self.assertNotIn('c.1G>G', df_scores[constants.nt_variant_col])
        self.assertNotIn('p.Ala1=', df_counts[constants.pro_variant_col])
        self.assertNotIn('p.Ala1=', df_scores[constants.pro_variant_col])


class TestEnrich2LoadInput(TestCase):
    def test_error_file_not_h5_or_tsv(self):
        path = os.path.join(DATA_DIR, 'empiric.xlsx')
        p = enrich2.Enrich2(path, wt_sequence='AAA')
        with self.assertRaises(TypeError):
            p.load_input_file()
            
    def test_scores_tsv_missing_score_column(self):
        path = os.path.join(DATA_DIR, 'enrich2.tsv')
        p = enrich2.Enrich2(
            path, wt_sequence='AAA',
            score_column='scores',
            hgvs_column='sequence',
            input_type=constants.score_type
        )
        with self.assertRaises(KeyError):
            p.load_input_file()

    def test_input_type_counts_doesnt_raise_keyerror(self):
        path = os.path.join(DATA_DIR, 'enrich2.tsv')
        p = enrich2.Enrich2(
            path, wt_sequence='AAA',
            hgvs_column='sequence',
            input_type=constants.count_type
        )
        p.load_input_file()

    def test_scores_tsv_missing_hgvs_column(self):
        path = os.path.join(DATA_DIR, 'enrich2.tsv')
        p = enrich2.Enrich2(
            path, wt_sequence='AAA',
            hgvs_column='hgvs'
        )
        with self.assertRaises(KeyError):
            p.load_input_file()


class TestEnrich2ParseRow(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(DATA_DIR, 'dummy.h5')
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence='ACT')
        self.bin.append(self.path.replace('.h5', ''))

    def test_valueerr_cannot_deduce_type(self):
        with self.assertRaises(exceptions.InvalidVariantType):
            self.enrich2.parse_row(('c.1_2del', None))

    def test_delegates_to_dna(self):
        for prefix in hgvsp.constants.dna_prefix:
            variant = '{0}.1A>G,{0}.2C>G'.format(prefix)
            expected = '{0}.[1A>G;2C>G]'.format(prefix), None
            self.assertEqual(expected, self.enrich2.parse_row((variant, None)))

    def test_delegates_to_protein(self):
        variant = '{0}.Thr1=,{0}.Thr1Gly'.format(
            hgvsp.constants.protein_prefix)
        expected = None, '{0}.[Thr1=;Thr1Gly]'.format(
            hgvsp.constants.protein_prefix)
        self.assertEqual(expected, self.enrich2.parse_row((variant, None)))

    def test_nt_variant_is_none_special_variant_is_from_synonymous_table(self):
        self.assertEqual(
            (None, constants.enrich2_synonymous),
            self.enrich2.parse_row(
                (constants.enrich2_synonymous, constants.synonymous_table)
            ))
    
    @mock.patch("mavedbconvert.enrich2.apply_offset",
                return_value='c.3T>C (p.Thr1=)')
    def test_calls_apply_offset_to_variant(self, patch):
        variant = 'c.3T>C (p.=)'
        self.enrich2.parse_row((variant, None))
        patch.assert_called()
        
    def test_delegate_to_multi(self):
        variant = 'c.3T>C (p.Thr1=)'
        expected = ('c.3T>C', 'p.Thr1=')
        self.assertEqual(expected, self.enrich2.parse_row((variant, None)))

    def test_returns_special_variant_as_tuple_non_synonymous_table(self):
        self.assertEqual(self.enrich2.parse_row(('_wt', None)), ('_wt', '_wt'))
        self.assertEqual(self.enrich2.parse_row(('_sy', None)), ('_sy', '_sy'))

    def test_strips_whitespace(self):
        self.assertEqual(
            self.enrich2.parse_row((' c.1A>G ', None)), ('c.1A>G', None))
        
    def test_converts_X_to_N(self):
        self.assertEqual(
            self.enrich2.parse_row(('c.1A>X', None)), ('c.1A>N', None))
        
    def test_converts_triple_q_to_single_q(self):
        self.assertEqual(
            self.enrich2.parse_row(('p.Thr1???', None)), (None, 'p.Thr1?'))


# Protein parsing tests
# --------------------------------------------------------------------------- #
class TestProteinHGVSParsing(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(DATA_DIR, 'dummy.h5')
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence='AAA')
        self.bin.append(self.path.replace('.h5', ''))

    def test_combines_into_multi_variant_syntax(self):
        result = self.enrich2.parse_protein_variant('p.L5G, p.L6G')
        self.assertEqual(result, 'p.[L5G;L6G]')

    def test_combines_list_into_multi_variant_syntax(self):
        result = self.enrich2.parse_protein_variant(['p.L5G', 'p.L6G'])
        self.assertEqual(result, 'p.[L5G;L6G]')

    def test_does_not_squash_single_variant(self):
        result = self.enrich2.parse_protein_variant('p.L5G')
        self.assertEqual(result, 'p.L5G')

    def test_passes_on_sy_or_wt(self):
        self.assertEqual(self.enrich2.parse_protein_variant('_wt'), '_wt')
        self.assertEqual(self.enrich2.parse_protein_variant('_sy'), '_sy')

    def test_valueerr_not_a_valid_protein_syntax(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant('c.101A>G')
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant('p.101A>G')
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant('random')
        with self.assertRaises(ValueError):
            self.enrich2.parse_protein_variant('(random)')

    def test_removes_brackets(self):
        result = self.enrich2.parse_protein_variant('(p.L4G),(p.L5G)')
        self.assertEqual(result, 'p.[L4G;L5G]')

        result = self.enrich2.parse_protein_variant('(p.L4G)')
        self.assertEqual(result, 'p.L4G')

        result = self.enrich2.parse_protein_variant(['(p.L4G)'])
        self.assertEqual(result, 'p.L4G')

    def test_strips_ws(self):
        result = self.enrich2.parse_protein_variant(' p.Gly5Leu ')
        self.assertEqual(result, 'p.Gly5Leu')
        result = self.enrich2.parse_protein_variant(' p.L4G, p.L5G ')
        self.assertEqual(result, 'p.[L4G;L5G]')

    def test_removes_duplicates(self):
        result = self.enrich2.parse_protein_variant('p.L5G, p.L5G')
        self.assertEqual(result, 'p.L5G')

    def test_maintains_ordering(self):
        result = self.enrich2.parse_protein_variant('p.L5G, p.L4G')
        self.assertEqual(result, 'p.[L5G;L4G]')


# Nucleotide parsing tests
# --------------------------------------------------------------------------- #
class TestNucleotideHGVSParing(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(DATA_DIR, 'dummy.h5')
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence='AAA')
        self.bin.append(self.path.replace('.h5', ''))
        
    def test_parses_non_coding_nt_variants_into_multi_variant(self):
        nt = self.enrich2.parse_nucleotide_variant(
            'n.-455T>A, n.-122A>T, n.-101A>T, n.-42T>A')
        self.assertEqual(nt, 'n.[-455T>A;-122A>T;-101A>T;-42T>A]')

    def test_combines_into_multi_variant_syntax(self):
        result = self.enrich2.parse_nucleotide_variant('c.2A>G,c.1A>G')
        self.assertEqual(result, 'c.[2A>G;1A>G]')

    def test_combines_list_into_multi_variant_syntax(self):
        result = self.enrich2.parse_nucleotide_variant(['c.2A>G', 'c.1A>G'])
        self.assertEqual(result, 'c.[2A>G;1A>G]')

    def test_does_not_squash_single_variant(self):
        result = self.enrich2.parse_nucleotide_variant('c.2A>G')
        self.assertEqual(result, 'c.2A>G')

    def test_passes_on_sy_or_wt(self):
        self.assertEqual(self.enrich2.parse_nucleotide_variant('_wt'), '_wt')
        self.assertEqual(self.enrich2.parse_nucleotide_variant('_sy'), '_sy')

    def test_valueerr_not_a_valid_protein_syntax(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant('p.101A>G')
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant('p.L5G')
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant('random')
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant('()')
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant('c.[2A>G;1A>G]')

    def test_valueerr_multi_prefix_types(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_nucleotide_variant('c.1A>G;n.2A>G')

    def test_strips_ws(self):
        result = self.enrich2.parse_nucleotide_variant(' c.101A>G ')
        self.assertEqual(result, 'c.101A>G')
        result = self.enrich2.parse_nucleotide_variant(' c.2A>G, c.2A>G ')
        self.assertEqual(result, 'c.[2A>G;2A>G]')


# Mixed parsing tests
# --------------------------------------------------------------------------- #
class TestEnrich2MixedHGVSParsing(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(DATA_DIR, 'dummy.h5')
        self.wt = 'ACT'
        self.wt_aa = constants.AA_CODES[constants.CODON_TABLE[self.wt]]
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence=self.wt)
        self.bin.append(self.path.replace('.h5', ''))

    def test_parses_nt_variants_into_multi_variant(self):
        nt, _ = self.enrich2.parse_mixed_variant(
            'c.1A>T (p.Thr1Tyr), c.2C>A (p.Thr1Tyr)',)
        self.assertEqual(nt, 'c.[1A>T;2C>A]')
        self.assertIsNotNone(hgvsp.multi_variant_re.fullmatch(nt))

    def test_parses_pro_variants_into_multi_variant(self):
        self.enrich2.wt_sequence = 'ACTCAA'
        _, pro = self.enrich2.parse_mixed_variant(
            'c.1A>T (p.Thr1Pro), c.4C>A (p.Gln2Lys)')
        self.assertEqual(pro, 'p.[Thr1Pro;Gln2Lys]')
        self.assertIsNotNone(hgvsp.multi_variant_re.fullmatch(pro))

    def test_error_list_len_different(self):
        # ValueError when attempting tuple unpacking
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant('c.1A>G, c.2T>A (p.Lys1Arg)')
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant('p.Lys4Arg, c.2T>A (p.Lys1Arg)')

    def test_variant_order_maintained(self):
        self.enrich2.wt_sequence = 'AAAAAT'
        nt, pro = self.enrich2.parse_mixed_variant(
            'c.1= (p.Lys1Ile), c.6T>G (p.Asn2Lys), c.2A>T (p.Lys1Ile)')
        self.assertEqual(nt, 'c.[1=;6T>G;2A>T]')
        self.assertEqual(pro, 'p.[Lys1Ile;Asn2Lys]')

    @mock.patch.object(enrich2.Enrich2, 'infer_silent_aa_substitution',
                       return_value='p.Lys1=')
    def test_groups_codons(self, patch):
        self.enrich2.wt_sequence = 'AAAAAT'
        variant = 'c.1= (p.=), c.6T>G (p.Asn2Lys), c.2= (p.=)'
        _, _ = self.enrich2.parse_mixed_variant(variant)
        patch.assert_called_with(*(['c.1=', 'c.2='], variant))

    @mock.patch.object(enrich2.Enrich2, 'infer_silent_aa_substitution',
                       return_value='p.Lys1=')
    def test_calls_infer_with_synonymous_variants_only(self, patch):
        self.enrich2.wt_sequence = 'AAAAAT'
        variant = 'c.1= (p.=), c.6T>G (p.Asn2Lys), c.2= (p.Lys1=)'
        _, _ = self.enrich2.parse_mixed_variant(variant)
        patch.assert_called_with(*(['c.1=', ], variant))

    def test_nt_variant_is_none_special_variant_is_from_synonymous_table(self):
        self.assertEqual(
            (None, constants.enrich2_synonymous),
            self.enrich2.parse_row(
                (constants.enrich2_synonymous, constants.synonymous_table)
            ))

    def test_valueerror_multiple_prefix_types(self):
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant(
                'c.1A>G (p.=), r.2u>a (p.Lys4Arg)')
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant(
                'c.1A>G (p.=), n.2T>A (p.Lys4Arg)')
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant(
                'c.1A>G (p.=), (p.Lys4Arg)')
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant(
                'c.1A>G (p.=), p.Lys4Arg')
        with self.assertRaises(ValueError):
            self.enrich2.parse_mixed_variant(
                'c.1A>G (p.=), c.2T>A (g.Lys4Arg)')

    def test_doesnt_collapse_single_variants_into_multivariant(self):
        nt, pro = self.enrich2.parse_mixed_variant('c.3T>C (p.=)')
        self.assertEqual(nt, 'c.3T>C')
        self.assertEqual(pro, 'p.Thr1=')
        self.assertIsNotNone(hgvsp.single_variant_re.fullmatch(nt))
        self.assertIsNotNone(hgvsp.single_variant_re.fullmatch(pro))

    def test_protein_set_as_nt_when_table_is_not_syn_and_variant_is_special(self):
        nt, pro = self.enrich2.parse_mixed_variant('_wt')
        self.assertEqual(nt, '_wt')
        self.assertEqual(pro, '_wt')

        nt, pro = self.enrich2.parse_mixed_variant('_sy')
        self.assertEqual(nt, '_sy')
        self.assertEqual(pro, '_sy')


class TestInferSilentAASub(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(DATA_DIR, 'dummy.h5')
        self.enrich2 = enrich2.Enrich2(self.path, wt_sequence='AAA', offset=0)
        self.bin.append(self.path.replace('.h5', ''))

    def test_valueerror_not_a_sub_event(self):
        with self.assertRaises(exceptions.InvalidVariantType):
            self.enrich2.infer_silent_aa_substitution('c.100_102del')

    def test_index_error_variant_pos_is_out_of_bounds_relative_to_wt(self):
        with self.assertRaises(IndexError):
            self.enrich2.infer_silent_aa_substitution('c.100A>G')

    def test_valueerror_base_in_wt_does_not_match_base_in_hgvs(self):
        self.enrich2.wt_sequence = 'TGA'
        with self.assertRaises(ValueError):
            self.enrich2.infer_silent_aa_substitution('c.1A>G')

    def test_correct_wt_aa_inferred(self):
        self.enrich2.wt_sequence = 'TCT'
        self.assertEqual(
            'p.Ser1=', self.enrich2.infer_silent_aa_substitution('c.3T>C')
        )

    def test_correct_aa_position_inferred(self):
        self.enrich2.wt_sequence = 'AAAGGGTCT'
        self.assertEqual(
            'p.Ser3=',
            self.enrich2.infer_silent_aa_substitution('c.9T>C')
        )

    def test_error_mutant_codon_does_not_match_wild_type(self):
        self.enrich2.wt_sequence = 'ATG'
        with self.assertRaises(ValueError):
            self.enrich2.infer_silent_aa_substitution('c.1A>C')

    def test_correctly_infers_aa_from_codon_group(self):
        self.enrich2.wt_sequence = 'TTA'
        group = ['c.1T>C', 'c.2=', 'c.3A>T']
        self.assertEqual(
            'p.Leu1=', self.enrich2.infer_silent_aa_substitution(group)
        )

    def test_valueerror_mixed_codons_in_group(self):
        with self.assertRaises(ValueError):
            self.enrich2.infer_silent_aa_substitution(['c.1T>C', 'c.5T>C'])

    def test_correctly_infers_aa_from_silent_variants(self):
        self.enrich2.wt_sequence = 'TTA'
        group = ['c.1=', 'c.2=', 'c.3=']
        self.assertEqual(
            'p.Leu1=', self.enrich2.infer_silent_aa_substitution(group)
        )


class TestApplyOffset(TestCase):
    def test_mixed_variant_uses_nt_position_to_compute_codon_pos(self):
        variant = 'c.-9A>T (p.Thr2Pro), c.-6C>A (p.Gln3Lys)'
        offset = -10
        self.assertEquals(
            'c.1A>T (p.Thr1Pro), c.4C>A (p.Gln2Lys)',
            enrich2.apply_offset(variant, offset)
        )

    def test_error_position_after_offset_non_positive(self):
        with self.assertRaises(ValueError):
            enrich2.apply_offset('c.1A>T', 10)
            
        with self.assertRaises(ValueError):
            enrich2.apply_offset('p.Leu1=', 10)
        
    def test_applies_offset_to_non_mixed_variant(self):
        variant = 'n.-455T>A, n.-122A>T, n.-101A>T, n.-42T>A'
        offset = -456
        self.assertEquals(
            'n.1T>A, n.334A>T, n.355A>T, n.414T>A',
            enrich2.apply_offset(variant, offset)
        )
        self.assertEquals(
            'n.1T>A',
            enrich2.apply_offset('n.-455T>A', offset)
        )
        
    def test_applies_offset_to_protein_variant_modulo_3(self):
        variant = 'p.Leu10=, p.Leu13='
        offset = 10
        self.assertEquals(
            'p.Leu7=, p.Leu10=', enrich2.apply_offset(variant, offset)
        )
        self.assertEquals(
            'p.Leu7=', enrich2.apply_offset('p.Leu10=', offset)
        )
    
    @mock.patch.object(enrich2.base.BaseProgram, 'validate_against_wt_sequence')
    def test_validates_against_wt_sequence(self, patch):
        variant = 'c.-9C>T'
        path = os.path.join(DATA_DIR, 'enrich1.tsv')
        p = enrich2.Enrich2(path, wt_sequence='ACT')
        enrich2.apply_offset(variant, offset=-10, enrich2=p)  # pass
        patch.assert_called_with(*('c.1C>T',))

    def test_value_error_base_mismatch_after_offset_applied(self):
        variant = 'c.-9G>T'
        path = os.path.join(DATA_DIR, 'enrich1.tsv')
        p = enrich2.Enrich2(path, wt_sequence='ACT')
        with self.assertRaises(ValueError):
            enrich2.apply_offset(variant, offset=-10, enrich2=p)
            
    @mock.patch.object(enrich2.base.BaseProgram, 'validate_against_protein_sequence')
    def test_validates_against_pro_sequence(self, patch):
        variant = 'p.Gly3Leu'
        path = os.path.join(DATA_DIR, 'enrich1.tsv')
        p = enrich2.Enrich2(path, wt_sequence='ACG')
        enrich2.apply_offset(variant, offset=6, enrich2=p)  # pass
        patch.assert_called_with(*('p.Gly1Leu',))

    def test_value_error_pro_mismatch_after_offset_applied(self):
        variant = 'p.Gly3Leu'
        path = os.path.join(DATA_DIR, 'enrich1.tsv')
        p = enrich2.Enrich2(path, wt_sequence='ACG')
        with self.assertRaises(ValueError):
            enrich2.apply_offset(variant, offset=6, enrich2=p)
