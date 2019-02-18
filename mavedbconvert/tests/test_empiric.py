import os
from unittest import TestCase

import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal, assert_series_equal

from .. import empiric, constants

from . import ProgramTestCase


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(BASE_DIR + '/data/')


class TestInferProEvent(TestCase):
    def test_infers_equal_event(self):
        self.assertEqual(
            empiric.infer_pro_substitution(
                mut_aa='V', wt_aa='v', codon_pos=0
            ), 'p.Val1='
        )
        
    def test_infers_sub_event_event(self):
        self.assertEqual(
            empiric.infer_pro_substitution(
                mut_aa='V', wt_aa='F', codon_pos=0
            ), 'p.Phe1Val'
        )
        
    def test_converts_triple_q_to_single_q(self):
        self.assertEqual(
            empiric.infer_pro_substitution(
                mut_aa='???', wt_aa='v', codon_pos=0
            ), 'p.Val1?'
        )


class TestInferNTEvent(TestCase):
    def test_infers_equal_event(self):
        self.assertEqual(
            empiric.infer_nt_substitution(
                wt_codon='aaa', mut_codon='AAA', codon_pos=0
            ), 'c.[1=;2=;3=]'
        )
    
    def test_infers_sub_event_event(self):
        self.assertEqual(
            empiric.infer_nt_substitution(
                wt_codon='ATC', mut_codon='GTA', codon_pos=0
            ), 'c.[1A>G;2=;3C>A]'
        )
        
    def test_adds_codon_pos_multiplied_by_3_to_position(self):
        self.assertEqual(
            empiric.infer_nt_substitution(
                wt_codon='ATC', mut_codon='GTA', codon_pos=1
            ), 'c.[4A>G;5=;6C>A]'
        )
        
    def test_converts_X_to_N(self):
        self.assertEqual(
            empiric.infer_nt_substitution(
                wt_codon='ATC', mut_codon='GTX', codon_pos=1
            ), 'c.[4A>G;5=;6C>N]'
        )


class TestEmpiric(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.input = os.path.join(DATA_DIR, 'empiric.xlsx')
        self.empiric = empiric.Empiric(
            src=self.input, wt_sequence='AAA', one_based=False,
        )

    def test_error_missing_amino_acid(self):
        for nan in constants.extra_na:
            df = pd.DataFrame({
                'Position': [0, ], 'Amino Acid': [nan, ], 'row_num': [0]
            })
            self.empiric.validate_columns(df)
            with self.assertRaises(ValueError):
                self.empiric.parse_row(row=df.iloc[0, :],)

    def test_value_error_codon_doesnt_match_aa_column(self):
        with self.assertRaises(ValueError):
            df = pd.DataFrame({
                'Position': [0, ], 'Amino Acid': ['V', ],
                'Codon': ['AAT', ], 'row_num': [0]
            })
            self.empiric.validate_columns(df)
            self.empiric.parse_row(row=df.iloc[0, :],)
                
    def test_error_infer_nt_true_but_missing_codon_value(self):
        for nan in constants.extra_na:
            df = pd.DataFrame({
                'Position': [0, ], 'Amino Acid': ['N', ],
                'row_num': [0], 'Codon': [nan, ],
            })
            self.empiric.validate_columns(df)
            with self.assertRaises(ValueError):
                self.empiric.parse_row(row=df.iloc[0, :],)

    def test_index_error_negative_position(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['K', ],
            'row_num': [0], 'Codon': ['AAA', ],
        })
        self.empiric.validate_columns(df)
        self.empiric.one_based = True
        with self.assertRaises(IndexError):
            self.empiric.parse_row(row=df.iloc[0, :],)

    def test_index_error_out_of_codon_bounds(self):
        df = pd.DataFrame({
            'Position': [56, ], 'Amino Acid': ['K', ],
            'row_num': [0], 'Codon': ['AAA', ],
        })
        self.empiric.validate_columns(df)
        with self.assertRaises(IndexError):
            self.empiric.parse_row(row=df.iloc[0, :], )
    
    def test_amino_acid_column_is_case_insensitive(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['v', ],
            'row_num': [0], 'Codon': ['GTA', ],
        })
        self.empiric.validate_columns(df)
        _, hgvs_pro = self.empiric.parse_row(row=df.iloc[0, :],)
        self.assertEqual(hgvs_pro, 'p.Lys1Val')
        
    def test_infers_hgvs_pro_event_from_one_based_position(self):
        df = pd.DataFrame({
            'Position': [1, ], 'Amino Acid': ['V', ],
            'row_num': [0], 'Codon': ['GTA', ],
        })
        self.empiric.validate_columns(df)
        self.empiric.one_based = True
        _, hgvs_pro = self.empiric.parse_row(row=df.iloc[0, :],)
        self.assertEqual(hgvs_pro, 'p.Lys1Val')

    def test_infers_hgvs_pro_event_from_zero_based_position(self):
        df = pd.DataFrame({
            'Position': [1, ], 'Amino Acid': ['V', ],
            'row_num': [0], 'Codon': ['GTA', ],
        })
        self.empiric.validate_columns(df)
        self.empiric.wt_sequence = 'GTAAAA'
        self.empiric.one_based = False
        _, hgvs_pro = self.empiric.parse_row(row=df.iloc[0, :],)
        self.assertEqual(hgvs_pro, 'p.Lys2Val')

    def test_protein_output_is_singular_when_inferring_nt(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['V', ],
            'row_num': [0], 'Codon': ['GTA', ],
        })
        self.empiric.validate_columns(df)
        hgvs_nt, hgvs_pro = self.empiric.parse_row(row=df.iloc[0, :],)
        self.assertEqual(hgvs_nt, 'c.[1A>G;2A>T;3=]')
        self.assertEqual(hgvs_pro, 'p.Lys1Val')
        
    def test_hgvs_nt_is_none_when_codon_is_not_in_axes(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['V', ],
            'row_num': [0],
        })
        self.empiric.validate_columns(df)
        hgvs_nt, _ = self.empiric.parse_row(row=df.iloc[0, :],)
        self.assertIsNone(hgvs_nt)

    def test_correctly_infers_hgvs_nt_positions_when_zero_based(self):
        df = pd.DataFrame({
            'Position': [1, ], 'Amino Acid': ['V', ],
            'row_num': [0], 'Codon': ['GTA', ],
        })
        self.empiric.validate_columns(df)
        self.empiric.one_based = False
        self.empiric.wt_sequence = 'GGGAAT'
        hgvs_nt, _ = self.empiric.parse_row(row=df.iloc[0, :],)
        self.assertEqual(hgvs_nt, 'c.[4A>G;5A>T;6T>A]')

    def test_correctly_infers_hgvs_nt_positions_when_one_based(self):
        df = pd.DataFrame({
            'Position': [1, ], 'Amino Acid': ['N', ],
            'Codon': ['AAT', ],
        })
        self.empiric.validate_columns(df)
        self.empiric.one_based = True
        self.empiric.wt_sequence = 'GTAAAA'
        hgvs_nt, _ = self.empiric.parse_row(row=df.iloc[0, :], )
        self.assertEqual(hgvs_nt, 'c.[1G>A;2T>A;3A>T]')


class TestEmpiricValidateColumns(TestCase):
    def setUp(self):
        super().setUp()
        self.input = os.path.join(DATA_DIR, 'empiric.xlsx')
        self.empiric = empiric.Empiric(
            src=self.input, wt_sequence='AAA',
            one_based=False,
        )
        
    def test_error_cannot_find_case_insensitive_aa_column(self):
        df = pd.DataFrame({
            'Position': [1, ], 'aa': ['N', ], 'Codon': ['AAT', ],
        })
        with self.assertRaises(ValueError):
            self.empiric.validate_columns(df)

    def test_error_cannot_find_case_insensitive_position_column(self):
        df = pd.DataFrame({
            'pos': [1, ], 'Amino Acid': ['N', ], 'Codon': ['AAT', ],
        })
        with self.assertRaises(ValueError):
            self.empiric.validate_columns(df)
            
    def test_sets_codon_column_as_none_if_not_present(self):
        df = pd.DataFrame({
            'Position': [1, ], 'Amino Acid': ['N', ],
        })
        self.empiric.validate_columns(df)
        self.assertEqual(self.empiric.codon_column, None)
        
    def test_sets_codon_column_if_present(self):
        df = pd.DataFrame({
            'Position': [1, ], 'Amino Acid': ['N', ], 'Codon': ['AAT', ],
        })
        self.empiric.validate_columns(df)
        self.assertEqual(self.empiric.codon_column, "Codon")
        
    def test_sets_position_column(self):
        df = pd.DataFrame({
            'Position': [1, ], 'Amino Acid': ['N', ], 'Codon': ['AAT', ],
        })
        self.empiric.validate_columns(df)
        self.assertEqual(self.empiric.position_column, "Position")
        
    def test_sets_aa_column(self):
        df = pd.DataFrame({
            'Position': [1, ], 'amino acid': ['N', ], 'Codon': ['AAT', ],
        })
        self.empiric.validate_columns(df)
        self.assertEqual(self.empiric.aa_column, "amino acid")


class TestEmpiricParseInput(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.input = os.path.join(DATA_DIR, 'empiric.xlsx')
        self.empiric = empiric.Empiric(
            src=self.input, wt_sequence='AAA',
            one_based=False, input_type='scores',
            score_column='A',
        )

    def test_deletes_position_amino_acid_codon_row_num_columns(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['N', ], 'Codon': ['AAT', ],
            'A': [1.2, ]
        })
        result = self.empiric.parse_input(df)
        self.assertNotIn('Position', result.columns)
        self.assertNotIn('Amino Acid', result.columns)
        self.assertNotIn('Codon', result.columns)
        self.assertNotIn('row_num', result.columns)
        
    def test_keeps_additional_non_score_columns(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['N', ],
            'Codon': ['AAT', ], 'A': [1.2], 'B': [2.4],
        })
        result = self.empiric.parse_input(df)
        self.assertIn('B', result.columns)

    def test_renames_score_column_to_score_and_drops_original(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['N', ],
            'Codon': ['AAT', ], 'A': [1.2], 'B': [2.4],
        })
        result = self.empiric.parse_input(df)
        self.assertListEqual(list(df['A']), list(result['score']))
        self.assertIn('B', result.columns)
        self.assertNotIn('A', result.columns)
        
    def test_sets_hgvs_pro_column(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['N', ],
            'Codon': ['AAT', ], 'A': [1.2], 'B': [2.4],
        })
        result = self.empiric.parse_input(df)
        self.assertEqual(
            result[constants.pro_variant_col].values[0], 'p.Lys1Asn')

    def test_correctly_infers_hgvs_nt_column_when_codon_column_present(self):
        df = pd.DataFrame({
            'Position': [1, ], 'Amino Acid': ['N', ],
            'Codon': ['AAT', ], 'A': [1.2], 'B': [2.4],
        })
        self.empiric.one_based = False
        self.empiric.wt_sequence = 'GGGAAA'
        result = self.empiric.parse_input(df)
        self.assertEqual(
            result[constants.nt_variant_col].values[0], 'c.[4=;5=;6A>T]')

    def test_orders_columns(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['N', ],
            'Codon': ['AAT', ], 'A': [1.2], 'B': [2.4],
        })
        result = self.empiric.parse_input(df)
        self.assertEqual(list(result.columns).index(
            constants.nt_variant_col), 0)
        self.assertEqual(list(result.columns).index(
            constants.pro_variant_col), 1)

    def test_removes_null_columns(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['N', ],
            'Codon': ['AAT', ], 'B': [None,], 'A': [2.4],
        })
        result = self.empiric.parse_input(df)
        self.assertNotIn('B', result.columns)
        
    def test_drops_nt_when_codon_column_is_not_provided(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['N', ],
            'A': [1.2,], 'B': [2.4],
        })
        result = self.empiric.parse_input(df)
        self.assertNotIn(constants.nt_variant_col, result.columns)

    def test_drops_non_numeric_columns(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['N', ],
            'Codon': ['AAT', ], 'A': [1.2, ], 'B': ['a', ],
        })
        result = self.empiric.parse_input(df)
        self.assertNotIn('B', result.columns)

    def test_keeps_int_type_as_int(self):
        df = pd.DataFrame({
            'Position': [0, ], 'Amino Acid': ['N', ],
            'Codon': ['AAT', ], 'A': [1, ],
        })
        result = self.empiric.parse_input(df)
        self.assertTrue(np.issubdtype(
            result[constants.mavedb_score_column].values[0], np.signedinteger))


class TestEmpiricLoadInput(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(DATA_DIR, 'empiric.xlsx')
        self.tmp_path = os.path.join(DATA_DIR, 'tmp.csv')
        self.tmp_excel_path = os.path.join(DATA_DIR, 'tmp.xlsx')
        self.bin.append(self.tmp_path)

    def test_extra_na_load_as_nan(self):
        for value in constants.extra_na:
            df = pd.read_excel(self.path)
            df['A'] = [value, ] * len(df)
            df.to_csv(self.tmp_path, index=False, sep='\t',)
            e = empiric.Empiric(
                src=self.tmp_path, wt_sequence='TTTTCTTATTGT', 
                score_column='col_A', input_type=constants.score_type,
                one_based=False,
            )
            result = e.load_input_file()
            expected = pd.Series([np.NaN]*len(df), index=df.index, name='A')
            assert_series_equal(result['A'], expected)

    def test_loads_first_sheet_by_default(self):
        data = [
            {'Position': [0], 'Amino Acid': ['K'], 'score': [1.2, ]},
            {'Position': [1], 'Amino Acid': ['G'], 'score': [1.4, ]},
        ]
        self.mock_multi_sheet_excel_file(self.tmp_excel_path, data)
        p = empiric.Empiric(
            src=self.tmp_excel_path, wt_sequence='TTTTCTTATTGT', 
            score_column='score', input_type=constants.score_type)
        df = p.load_input_file()
        expected = pd.DataFrame(data[0])
        assert_frame_equal(df, expected)

    def test_handles_csv(self):
        df = pd.read_excel(self.path)
        df.to_csv(self.tmp_path, index=False, sep='\t',)
        e = empiric.Empiric(
            src=self.tmp_path, wt_sequence='TTTTCTTATTGT', 
            score_column='col_A', input_type=constants.score_type,
            one_based=False,
        )
        result = e.load_input_file()
        assert_frame_equal(result, df)

    def test_error_position_not_in_columns(self):
        df = pd.read_excel(self.path)
        df = df.drop(columns=['Position'])
        df.to_csv(self.tmp_path, index=False, sep='\t')
        with self.assertRaises(ValueError):
            e = empiric.Empiric(
                src=self.tmp_path, wt_sequence='TTTTCTTATTGT', 
                score_column='col_A', input_type=constants.score_type,
                one_based=False,
            )
            e.load_input_file()

    def test_error_amino_acid_not_in_columns(self):
        df = pd.read_excel(self.path)
        df = df.drop(columns=['Amino Acid'])
        df.to_csv(self.tmp_path, index=False, sep='\t')
        with self.assertRaises(ValueError):
            e = empiric.Empiric(
                src=self.tmp_path, wt_sequence='TTTTCTTATTGT', 
                score_column='col_A', input_type=constants.score_type,
                one_based=False,
            )
            e.load_input_file()

    def test_not_scores_column_but_input_type_is_scores(self):
        with self.assertRaises(ValueError):
            empiric.Empiric(
                src=self.tmp_path, wt_sequence='TTTTCTTATTGT', 
                score_column=None, input_type=constants.score_type,
                one_based=False,
            )
            
    def test_adds_offset_to_position_column(self):
        e = empiric.Empiric(
            src=self.path, wt_sequence='TTTTCTTATTGT',
            score_column='col_A', input_type=constants.score_type,
            one_based=False, offset=10
        )
        result = e.load_input_file()
        print(e.position_column)
        self.assertListEqual(
            list(result[e.position_column]),
            [10, 11, 12]
        )


class TestEmpiricConvert(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.path = os.path.join(DATA_DIR, 'empiric.xlsx')
        self.expected = os.path.join(DATA_DIR, 'empiric_expected.csv')
        self.empiric = empiric.Empiric(
            src=self.path, wt_sequence='TTTTCTTATTGT', 
            score_column='col_A', input_type=constants.score_type,
            one_based=False,
        )
        self.bin.append(self.empiric.output_file)

    def test_saves_to_dst(self):
        self.empiric.convert()
        self.assertTrue(os.path.isfile(self.empiric.output_file))

    def test_integration(self):
        self.empiric = empiric.Empiric(
            src=self.path, wt_sequence='TCTTATTGT',
            score_column='col_A', input_type=constants.score_type,
            one_based=False,
        )
        self.empiric.convert()
        assert_frame_equal(
            pd.read_csv(self.empiric.output_file, delimiter=','),
            pd.read_csv(self.expected, delimiter=','),
        )
