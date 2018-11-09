from unittest import TestCase

import numpy as np

from .. import utilities


class TestUtilities(TestCase):
    def test_slicer_returns_chunks_of_size_n(self):
        self.assertEqual(
            list(utilities.slicer('aaabbbccc', 3)),
            ['aaa', 'bbb', 'ccc']
        )

    def test_slicer_returns_clips_if_cannot_chunk(self):
        self.assertEqual(
            list(utilities.slicer('aaaabbbbccc', 4)),
            ['aaaa', 'bbbb', 'ccc']
        )


class TestTranslateWTSequence(TestCase):
    def test_translate_wt_seq_no_offset(self):
        self.assertEqual(
            utilities.translate_dna('GTGGCGGAG', offset=0),
            'VAE'
        )

    def test_translate_wt_seq_with_offset(self):
        self.assertEqual(
            utilities.translate_dna('GTGGCGGAG', offset=3),
            'AE'
        )

    def test_translate_wt_error_not_multiple_of_three(self):
        with self.assertRaises(ValueError):
            utilities.translate_dna('GTGG')

    def test_error_offset_negative(self):
        with self.assertRaises(ValueError):
            utilities.translate_dna('GTGGCGGAG', offset=-3)


class TestDNAtoAAPosition(TestCase):
    def test_error_pos_negative(self):
        with self.assertRaises(ValueError):
            utilities.dna_to_aa_position(-1)

    def test_error_zero_pos_one_based(self):
        with self.assertRaises(ValueError):
            utilities.dna_to_aa_position(0, True)

    def test_convert_dna_to_aa_position_0_based(self):
        self.assertEqual(1, utilities.dna_to_aa_position(0))
        self.assertEqual(1, utilities.dna_to_aa_position(1))
        self.assertEqual(1, utilities.dna_to_aa_position(2))

        self.assertEqual(2, utilities.dna_to_aa_position(3))
        self.assertEqual(2, utilities.dna_to_aa_position(4))
        self.assertEqual(2, utilities.dna_to_aa_position(5))

    def test_convert_dna_to_aa_position_1_based(self):
        self.assertEqual(1, utilities.dna_to_aa_position(1, True))
        self.assertEqual(1, utilities.dna_to_aa_position(2, True))
        self.assertEqual(1, utilities.dna_to_aa_position(3, True))

        self.assertEqual(2, utilities.dna_to_aa_position(4, True))
        self.assertEqual(2, utilities.dna_to_aa_position(5, True))
        self.assertEqual(2, utilities.dna_to_aa_position(6, True))


class TestFormatColumn(TestCase):
    def test_replaces_null_with_nan(self):
        self.assertIs(utilities.format_column(['   '])[0], np.NaN)
        self.assertIs(utilities.format_column(['none'])[0], np.NaN)
        self.assertIs(utilities.format_column(['NAN'])[0], np.NaN)
        self.assertIs(utilities.format_column(['na'])[0], np.NaN)
        self.assertIs(utilities.format_column(['undefined'])[0], np.NaN)

    def test_ignores_nan(self):
        self.assertIs(utilities.format_column([np.NaN])[0], np.NaN)

    def test_ignores_non_null(self):
        self.assertIs(utilities.format_column([0.0])[0], 0.0)

    def test_typecasts_to_float_by_default(self):
        self.assertIsInstance(utilities.format_column([1])[0], float)

    def test_replaces_null_with_none_if_astype_is_not_int_or_float(self):
        self.assertIs(utilities.format_column(['none'], astype=str)[0], None)


class TestIsNumeric(TestCase):
    def test_true_for_float(self):
        self.assertTrue(utilities.is_numeric(float))

    def test_true_for_int(self):
        self.assertTrue(utilities.is_numeric(int))

    def test_true_for_np_float(self):
        self.assertTrue(utilities.is_numeric(np.float))

    def test_true_for_np_int(self):
        self.assertTrue(utilities.is_numeric(np.int))

    def test_false_for_str(self):
        self.assertFalse(utilities.is_numeric(str))

    def test_false_for_object(self):
        self.assertFalse(utilities.is_numeric(object))

    def test_false_for_np_object(self):
        self.assertFalse(utilities.is_numeric(np.object))


class TestNucleotideSubstitutionEvent(TestCase):
    def test_error_invalid_dna_substitution_syntax(self):
        with self.assertRaises(ValueError):
            utilities.NucleotideSubstitutionEvent('c.100_101delins')

    def test_strips_ws(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent(' c.1A>G ').variant,
            'c.1A>G'
        )

    def test_parses_position(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.1A>G').position, 1)

    def test_parses_ref_base(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.1A>G').ref, 'A')
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.1=').ref, None)

    def test_parses_alt_base(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.1A>G').alt, 'G')
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.1=').alt, None)

    def test_infers_silent(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.1=').silent, True)

    def test_parses_prefix(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.1A>G').prefix, 'c')

    def test_formats_event_string_correctly(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.1A>G').event,
            '1A>G'
        )
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.1=').event,
            '1='
        )

    def test_infers_codon_position(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.3A>G').
                codon_position(one_based=False),
            2
        )
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.3A>G').
                codon_position(one_based=True),
            1
        )

    def test_infers_within_frame_position(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.3A>G').
                codon_frame_position(one_based=False),
            1
        )
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent('c.3A>G').
                codon_frame_position(one_based=True),
            3
        )


class TestProteinSubstitutionEvent(TestCase):
    def test_error_invalid_dna_substitution_syntax(self):
        with self.assertRaises(ValueError):
            utilities.ProteinSubstitutionEvent('p.100_101delins')

    def test_strips_ws(self):
        self.assertEqual(
            utilities.ProteinSubstitutionEvent(' p.Gly2Leu ').variant,
            'p.Gly2Leu'
        )

    def test_parses_position(self):
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.Gly2Leu').position, 2)

    def test_parses_ref_three_letter_aa(self):
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.G2L').ref, 'Gly')
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.Gly2Leu').ref, 'Gly')

    def test_parses_alt_three_letter_aa(self):
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.G2L').alt, 'Leu')
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.Gly2Leu').alt, 'Leu')

    def test_sets_alt_as_ref_in_silent_variant(self):
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.Gly2=').alt, 'Gly')

    def test_infers_silent(self):
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.Gly2=').silent, True)

    def test_parses_prefix(self):
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.Gly2=').prefix, 'p')

    def test_formats_event_string_correctly(self):
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.Gly2=').event, 'Gly2=')
        self.assertEqual(
            utilities.ProteinSubstitutionEvent('p.Gly2Leu').event, 'Gly2Leu')