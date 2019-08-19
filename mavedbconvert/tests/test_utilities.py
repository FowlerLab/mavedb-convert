from unittest import TestCase

import numpy as np

from .. import utilities, constants, exceptions


class TestSlicer(TestCase):
    def test_slicer_returns_chunks_of_size_n(self):
        self.assertEqual(list(utilities.slicer("aaabbbccc", 3)), ["aaa", "bbb", "ccc"])

    def test_slicer_returns_clips_if_cannot_chunk(self):
        self.assertEqual(
            list(utilities.slicer("aaaabbbbccc", 4)), ["aaaa", "bbbb", "ccc"]
        )


class TestTranslateWTSequence(TestCase):
    def test_translate_wt_seq_no_offset(self):
        self.assertEqual(utilities.translate_dna("GTGGCGGAG", offset=0), "VAE")

    def test_translate_wt_seq_with_offset(self):
        self.assertEqual(utilities.translate_dna("GTGGCGGAG", offset=3), "AE")

    def test_translate_wt_error_not_multiple_of_three(self):
        with self.assertRaises(ValueError):
            utilities.translate_dna("GTGG")

    def test_error_offset_negative(self):
        with self.assertRaises(ValueError):
            utilities.translate_dna("GTGGCGGAG", offset=-3)


class TestIsNull(TestCase):
    def test_is_null_true_for_none_nan_and_na(self):
        for v in constants.extra_na:
            self.assertTrue(utilities.is_null(v))

    def test_test_is_null_true_blank(self):
        self.assertTrue(utilities.is_null(" "))

    def test_is_null_case_insensitive(self):
        self.assertTrue(utilities.is_null("NONE"))

    def test_is_null_false(self):
        self.assertFalse(utilities.is_null("1.2"))


class TestFormatColumn(TestCase):
    def test_replaces_null_with_nan(self):
        self.assertIs(utilities.format_column(["   "])[0], np.NaN)
        self.assertIs(utilities.format_column(["none"])[0], np.NaN)
        self.assertIs(utilities.format_column(["NAN"])[0], np.NaN)
        self.assertIs(utilities.format_column(["na"])[0], np.NaN)
        self.assertIs(utilities.format_column(["undefined"])[0], np.NaN)

    def test_ignores_nan(self):
        self.assertIs(utilities.format_column([np.NaN])[0], np.NaN)

    def test_ignores_non_null(self):
        self.assertIs(utilities.format_column([0.0])[0], 0.0)

    def test_typecasts_to_float_by_default(self):
        self.assertIsInstance(utilities.format_column([1])[0], float)

    def test_replaces_null_with_none_if_astype_is_not_int_or_float(self):
        self.assertIs(utilities.format_column(["none"], astype=str)[0], None)


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
    def test_parses_negative_positions(self):
        nt = utilities.NucleotideSubstitutionEvent("n.-100A>T")
        self.assertEqual(nt.position, -100)

    def test_negative_rna_position_error(self):
        with self.assertRaises(IndexError):
            utilities.NucleotideSubstitutionEvent("r.-100a>u")

    def test_error_invalid_dna_substitution_syntax(self):
        with self.assertRaises(exceptions.InvalidVariantType):
            utilities.NucleotideSubstitutionEvent("c.100_101delins")

    def test_strips_ws(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent(" c.1A>G ").variant, "c.1A>G"
        )

    def test_parses_position(self):
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.1A>G").position, 1)
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.-1A>G").position, -1)

    def test_error_negative_position_when_computing_codon_position(self):
        nt = utilities.NucleotideSubstitutionEvent("c.-1A>G")
        with self.assertRaises(ValueError):
            nt.codon_position()
        with self.assertRaises(ValueError):
            nt.codon_frame_position()

    def test_parses_ref_base(self):
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.1A>G").ref, "A")
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.1=").ref, None)

    def test_parses_alt_base(self):
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.1A>G").alt, "G")
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.1=").alt, None)

    def test_infers_silent(self):
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.1=").silent, True)

    def test_parses_prefix(self):
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.1A>G").prefix, "c")

    def test_formats_event_string_correctly(self):
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.1A>G").event, "1A>G")
        self.assertEqual(utilities.NucleotideSubstitutionEvent("c.1=").event, "1=")

    def test_infers_codon_position(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent("c.3A>G").codon_position(
                one_based=False
            ),
            2,
        )
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent("c.3A>G").codon_position(
                one_based=True
            ),
            1,
        )

    def test_infers_within_frame_position(self):
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent("c.3A>G").codon_frame_position(
                one_based=False
            ),
            1,
        )
        self.assertEqual(
            utilities.NucleotideSubstitutionEvent("c.3A>G").codon_frame_position(
                one_based=True
            ),
            3,
        )


class TestProteinSubstitutionEvent(TestCase):
    def test_error_set_position_less_than_1(self):
        pro = utilities.ProteinSubstitutionEvent("p.Gly4Leu")
        with self.assertRaises(ValueError):
            pro.position -= 4

    def test_error_invalid_dna_substitution_syntax(self):
        with self.assertRaises(exceptions.InvalidVariantType):
            utilities.ProteinSubstitutionEvent("p.100_101delins")

    def test_strips_ws(self):
        self.assertEqual(
            utilities.ProteinSubstitutionEvent(" p.Gly2Leu ").variant, "p.Gly2Leu"
        )

    def test_parses_position(self):
        self.assertEqual(utilities.ProteinSubstitutionEvent("p.Gly2Leu").position, 2)

    def test_parses_ref_three_letter_aa(self):
        self.assertEqual(utilities.ProteinSubstitutionEvent("p.G2L").ref, "Gly")
        self.assertEqual(utilities.ProteinSubstitutionEvent("p.Gly2Leu").ref, "Gly")

    def test_parses_alt_three_letter_aa(self):
        self.assertEqual(utilities.ProteinSubstitutionEvent("p.G2L").alt, "Leu")
        self.assertEqual(utilities.ProteinSubstitutionEvent("p.Gly2Leu").alt, "Leu")

    def test_sets_alt_as_ref_in_silent_variant(self):
        self.assertEqual(utilities.ProteinSubstitutionEvent("p.Gly2=").alt, "Gly")

    def test_infers_silent(self):
        self.assertEqual(utilities.ProteinSubstitutionEvent("p.Gly2=").silent, True)

    def test_parses_prefix(self):
        self.assertEqual(utilities.ProteinSubstitutionEvent("p.Gly2=").prefix, "p")

    def test_formats_event_string_correctly(self):
        self.assertEqual(utilities.ProteinSubstitutionEvent("p.Gly2=").event, "Gly2=")
        self.assertEqual(
            utilities.ProteinSubstitutionEvent("p.Gly2Leu").event, "Gly2Leu"
        )


class TestSplitVariant(TestCase):
    def test_split_hgvs_singular_list_non_multi_variant(self):
        self.assertListEqual(["c.100A>G"], utilities.split_variant("c.100A>G"))

    def test_split_hgvs_returns_list_of_single_variants(self):
        self.assertListEqual(
            ["c.100A>G", "c.101A>G"], utilities.split_variant("c.[100A>G;101A>G]")
        )


class TestNormalizeVariant(TestCase):
    def test_stripts_white_space(self):
        self.assertEqual(utilities.normalize_variant(" c.1A>G "), "c.1A>G")

    def test_passes_on_none(self):
        self.assertIsNone(utilities.normalize_variant(None))

    def test_passes_on_special(self):
        self.assertEqual(utilities.normalize_variant("_wt"), "_wt")
        self.assertEqual(utilities.normalize_variant("_sy"), "_sy")

    def test_replaces_qmark_with_Xaa_in_protein_variant(self):
        self.assertEqual(utilities.normalize_variant("p.G4???"), "p.G4Xaa")
        self.assertEqual(utilities.normalize_variant("p.G4?"), "p.G4X")

    def test_ignores_invalid_protein_variant(self):
        self.assertEqual(utilities.normalize_variant("p.G4??"), "p.G4??")

    def test_replaces_X_with_N_in_dna_variant(self):
        for p in "cgnm":
            self.assertEqual(
                utilities.normalize_variant("{}.100A>X".format(p)),
                "{}.100A>N".format(p),
            )

        for p in "cgnm":
            self.assertEqual(
                utilities.normalize_variant("{}.100_102delinsXXX".format(p)),
                "{}.100_102delinsNNN".format(p),
            )

    def test_replaces_X_with_N_in_rna_variant(self):
        self.assertEqual(utilities.normalize_variant("r.100a>x"), "r.100a>n")
        self.assertEqual(
            utilities.normalize_variant("r.100_102delinsnnn"), "r.100_102delinsnnn"
        )


class TestFormatVariant(TestCase):
    def test_stripts_white_space(self):
        self.assertEqual(utilities.format_variant(" c.1A>G "), "c.1A>G")

    def test_passes_on_none(self):
        self.assertIsNone(utilities.format_variant(None))


class TestHGVSProFromEventList(TestCase):
    def test_returns_single_event(self):
        result = utilities.hgvs_pro_from_event_list(["L4V"])
        self.assertEqual(result, "p.L4V")

    def test_strips_whitespace(self):
        result = utilities.hgvs_pro_from_event_list([" L4V "])
        self.assertEqual(result, "p.L4V")

    def test_combines_muilt_events(self):
        result = utilities.hgvs_pro_from_event_list(["L4V", "G5*"])
        self.assertEqual(result, "p.[L4V;G5*]")

    def test_removes_duplicate_events(self):
        result = utilities.hgvs_pro_from_event_list(["L4V", "L4V"])
        self.assertEqual(result, "p.L4V")
        result = utilities.hgvs_pro_from_event_list(["L4V", "L4V", "L5V"])
        self.assertEqual(result, "p.[L4V;L5V]")

    def test_retains_ordering(self):
        result = utilities.hgvs_pro_from_event_list(["L5V", "L4V"])
        self.assertEqual(result, "p.[L5V;L4V]")

    def test_error_invalid_hgvs(self):
        with self.assertRaises(exceptions.HGVSMatchError):
            utilities.hgvs_pro_from_event_list(["aaaa"])


class TestHGVSNTFromEventList(TestCase):
    def test_returns_single_event(self):
        result = utilities.hgvs_nt_from_event_list(["45A>G"], prefix="c")
        self.assertEqual(result, "c.45A>G")

    def test_strips_whitespace(self):
        result = utilities.hgvs_nt_from_event_list([" 45A>G "], prefix="c")
        self.assertEqual(result, "c.45A>G")

    def test_combines_muilt_events(self):
        result = utilities.hgvs_nt_from_event_list(
            ["45A>G", "127_128delinsAAA"], prefix="n"
        )
        self.assertEqual(result, "n.[45A>G;127_128delinsAAA]")

    def test_keeps_duplicate_events(self):
        result = utilities.hgvs_nt_from_event_list(["45a>u", "45a>u"], prefix="r")
        self.assertEqual(result, "r.[45a>u;45a>u]")

    def test_error_invalid_hgvs(self):
        with self.assertRaises(exceptions.HGVSMatchError):
            utilities.hgvs_nt_from_event_list(["aaaa"], prefix="c")


class TestNonHgvsColumns(TestCase):
    def test_returns_non_hgvs_columns(self):
        self.assertListEqual(
            ["score"],
            list(
                utilities.non_hgvs_columns(
                    ["score", constants.nt_variant_col, constants.pro_variant_col]
                )
            ),
        )


class TestHgvsColumns(TestCase):
    def test_returns_only_hgvs_columns(self):
        self.assertListEqual(
            [constants.nt_variant_col, constants.pro_variant_col],
            list(
                utilities.hgvs_columns(
                    ["score", constants.nt_variant_col, constants.pro_variant_col]
                )
            ),
        )
