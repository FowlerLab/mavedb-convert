import os
import mock

from .. import base

from . import ProgramTestCase


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(BASE_DIR + '/data/')


class TestBaseProgram(ProgramTestCase):
    """
    Test __init__ correctly sets up read and write directories,
    sequence information etc.
    """
    def setUp(self):
        super().setUp()
        self.src = os.path.join(DATA_DIR, 'enrich1.tsv')
        self.src_with_spaces = os.path.join(DATA_DIR, 'enrich   1.tsv')
        self.h5_src = os.path.join(DATA_DIR, 'dummy.h5')

    def tearDown(self):
        for path in self.bin:
            if os.path.exists(path) and os.path.isfile(path):
                os.remove(path)
            elif os.path.exists(path) and os.path.isdir(path):
                os.removedirs(path)

    def test_sets_directory_as_input_directory_if_dst_is_none(self):
        p = base.BaseProgram(src=self.src, dst=None, wt_sequence='AAA', )
        self.assertEqual(p.dst, DATA_DIR)

    def test_error_file_not_readable(self):
        with self.assertRaises(IOError):
            base.BaseProgram(src='a file', dst=None, wt_sequence='AAA', )

    def test_expands_user_and_norms_dst(self):
        p = base.BaseProgram(src=self.src, dst='~//user//', wt_sequence='AAA',)
        self.assertEqual(
            p.dst,
            os.path.join(os.path.expanduser('~'), 'user')
        )

    def test_dir_with_input_fname_appended_when_h5_and_dst_is_none(self):
        p = base.BaseProgram(src=self.h5_src, dst=None, wt_sequence='AAA', )
        self.assertEqual(p.dst, os.path.join(DATA_DIR, 'dummy'))
        self.bin.append(os.path.join(DATA_DIR, 'dummy'))

    def test_creates_directory_tree_if_it_doesnt_exist(self):
        output = os.path.join(DATA_DIR, 'outer_dir/inner_dir/')
        base.BaseProgram(src=self.h5_src, dst=output, wt_sequence='AAA', )
        self.assertTrue(os.path.isdir(output))
        self.bin.append(output)

    @mock.patch('os.access')
    def test_checks_read_permission(self, patch):
        p = base.BaseProgram(src=self.src, dst=None, wt_sequence='AAA', )
        self.assertEqual(patch.call_args_list[0][0], (p.src, os.R_OK))

    @mock.patch('os.access')
    def test_checks_write_permission(self, patch):
        p = base.BaseProgram(src=self.src, dst=None, wt_sequence='AAA', )
        self.assertEqual(patch.call_args_list[1][0], (p.dst, os.W_OK))

    def test_splits_src_into_filename_and_ext(self):
        p = base.BaseProgram(src=self.src, dst=None, wt_sequence='AAA', )
        self.assertEqual(p.src_filename, 'enrich1')
        self.assertEqual(p.ext, '.tsv')

    def test_lower_cases_ext(self):
        p = base.BaseProgram(
            src=self.src.replace('tsv', 'TSV'), wt_sequence='AAA', )
        self.assertEqual(p.ext, '.tsv')

    def test_value_error_negative_offset(self):
        with self.assertRaises(ValueError):
            base.BaseProgram(
                src=self.src, wt_sequence='ATCA', offset=-1)

    def test_dst_filename_replaces_whitespace_with_underscores(self):
        p = base.BaseProgram(src=self.src_with_spaces, wt_sequence='AAA')
        self.assertEqual(p.dst_filename, 'mavedb_enrich_1.csv')

    def test_output_file_joins_dst_and_dst_filename(self):
        p = base.BaseProgram(src=self.src, wt_sequence='AAA')
        self.assertEqual(
            p.output_file,
            os.path.join(DATA_DIR, 'mavedb_enrich1.csv')
        )

    def test_output_directory_expands_user_and_norms_path(self):
        p = base.BaseProgram(src=self.src, wt_sequence='AAA')
        p.dst = '~//user//'
        self.assertEqual(
            p.output_directory,
            os.path.join(os.path.expanduser('~'), 'user')
        )

    # --- Test property setters --- #
    def test_wt_setter_upper_cases_wt_sequence(self):
        p = base.BaseProgram(src=self.src, wt_sequence='AAA', )
        p.wt_sequence = 'ggg'
        self.assertEqual(p.wt_sequence, 'GGG')

    def test_wt_setter_clips_offset_bases_from_wt_sequence(self):
        p = base.BaseProgram(src=self.src, wt_sequence='AAA')
        p.offset = 3
        p.wt_sequence = 'ATGCGA'
        self.assertEqual(p.wt_sequence, 'CGA')

    def test_wt_setter_translates_wt_sequence_from_offset(self):
        p = base.BaseProgram(src=self.src, wt_sequence='AAA', )
        p.offset = 1
        p.wt_sequence = 'ATCA'
        self.assertEqual(p.protein_sequence, 'S')

    def test_wt_setter_creates_codons_from_wt_sequence_and_clips_offset(self):
        p = base.BaseProgram(src=self.src, wt_sequence='AAA', )
        p.offset = 1
        p.wt_sequence = 'ATCA'
        self.assertEqual(p.codons, ['TCA', ])

    def test_wt_setter_value_error_not_valid_wt_sequence(self):
        p = base.BaseProgram(src=self.src, wt_sequence='AAA')
        with self.assertRaises(ValueError):
            p.wt_sequence = 'fff'

    def test_wt_setter_can_set_wt_sequence_and_offset_with_tuple(self):
        p = base.BaseProgram(src=self.src, wt_sequence='AAA')
        p.wt_sequence = 'ATCA', 1
        self.assertEqual(p.offset, 1)
        self.assertEqual(p.wt_sequence, 'TCA')
        self.assertEqual(p.protein_sequence, 'S')
        self.assertEqual(p.codons, ['TCA', ])

    def test_value_error_offset_negative(self):
        p = base.BaseProgram(src=self.src, wt_sequence='AAA')
        with self.assertRaises(ValueError):
            p.offset = -1


class TestBaseProgramValidateAgainstWTSeq(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.src = os.path.join(DATA_DIR, 'enrich1.tsv')
        self.base = base.BaseProgram(
            src=self.src, wt_sequence='ATG', one_based=True)

    def test_error_not_a_dna_sub(self):
        with self.assertRaises(ValueError):
            self.base.validate_against_wt_sequence('p.Gly1Leu')

    def test_passes_on_special_and_silent(self):
        self.base.validate_against_wt_sequence('_wt')
        self.base.validate_against_wt_sequence('_sy')
        self.base.validate_against_wt_sequence('c.1=')

    def test_passes_when_reference_base_matches(self):
        self.base.one_based = False
        self.base.validate_against_wt_sequence('c.1T>G')

        self.base.one_based = True
        self.base.validate_against_wt_sequence('c.1A>G')

    def test_error_when_reference_base_doesnt_match_wt(self):
        with self.assertRaises(ValueError):
            self.base.one_based = False
            self.base.validate_against_wt_sequence('c.1C>G')

        with self.assertRaises(ValueError):
            self.base.one_based = True
            self.base.validate_against_wt_sequence('c.1T>G')

    def test_error_negative_position(self):
        with self.assertRaises(IndexError):
            self.base.one_based = True
            self.base.validate_against_wt_sequence('c.0T>G')

    def test_validates_multi(self):
        self.base.validate_against_wt_sequence('c.[1A>G;2T>G]')
        with self.assertRaises(ValueError):
            self.base.validate_against_wt_sequence('c.[1A>G;2A>G]')


class TestBaseProgramValidateAgainstProteinSeq(ProgramTestCase):
    def setUp(self):
        super().setUp()
        self.src = os.path.join(DATA_DIR, 'enrich1.tsv')
        self.base = base.BaseProgram(
            src=self.src, wt_sequence='ATGAAA', one_based=True)

    def test_error_not_a_protein_sub(self):
        with self.assertRaises(ValueError):
            self.base.validate_against_protein_sequence('c.1A>G')

    def test_passes_on_special_and_silent(self):
        self.base.validate_against_protein_sequence('_wt')
        self.base.validate_against_protein_sequence('_sy')
        self.base.validate_against_protein_sequence('p.=')

    def test_passes_when_reference_aa_matches(self):
        self.base.one_based = False
        self.base.validate_against_protein_sequence('p.Met0Lys')

        self.base.one_based = True
        self.base.validate_against_protein_sequence('p.Met1Lys')

    def test_error_when_reference_base_doesnt_match_wt(self):
        with self.assertRaises(ValueError):
            self.base.one_based = False
            self.base.validate_against_protein_sequence('p.Met1Lys')

        with self.assertRaises(ValueError):
            self.base.one_based = True
            self.base.validate_against_protein_sequence('p.Met2Lys')

    def test_error_negative_position(self):
        with self.assertRaises(IndexError):
            self.base.one_based = True
            self.base.validate_against_protein_sequence('p.Met0Lys')

    def test_validates_multi(self):
        self.base.validate_against_protein_sequence('p.[Met1Lys;Lys2=]')
        with self.assertRaises(ValueError):
            self.base.validate_against_protein_sequence('p.[Met1Lys;Met2=]')
