import os
import re
import logging
import numpy as np

from hgvsp import is_multi

from . import LOGGER, utilities, constants


logger = logging.getLogger(LOGGER)


__all__ = ['BaseProgram', ]


class BaseProgram(object):
    """
    Convert an input file to MaveDB_ compliant counts or scores files.

    Attributes
    ----------
    src : str
        Source file to convert.

    wt_sequence : str
        A DNA wild-type sequence used for validation and inference of
        missing variants.

    offset : int, optional.
        The number of bases to clip in `wt_seq`. The resulting sequence
        should be the coding sequence analyzed in your input file.

    dst : str, optional.
        Directory to save the output to. Inferred as input directory if not
        specified.

    one_based : bool, optional.
        Set to `True` if input positions are `1-based` relative to the
        the wild-type sequence specified.

    skip_header_rows : int
        The number of lines to skip at the start of the input file. Only
        applicable to `Excel` and `TSV` files.

    skip_header_rows : int
        The number of lines to skip at the end of the input file. Only
        applicable to `Excel` and `TSV` files.

    sheet_name : str, optional.
        Name of the sheet to convert in an excel file.

    score_column : str, optional.
        The name of the column in the input column which you would like to
        set as the MaveDB score column. Ignored if `input_type` is 'counts'.

    input_type : str, optional.
        The MaveDB file type. Can be either 'scores' or 'counts'.
    """
    def __init__(self, src, wt_sequence, offset=0, dst=None, one_based=True,
                 skip_header_rows=0, skip_footer_rows=0, score_column=None,
                 input_type=None, sheet_name=None):
        # Check the input is a readable file.
        self.src = os.path.normpath(os.path.expanduser(src))
        logger.info("Checking read permission for '{}'".format(self.src))
        os.access(self.src, os.R_OK)

        src_filename, ext = os.path.splitext(os.path.split(src)[1])
        self.src_filename = src_filename
        self.dst_filename = 'mavedb_{}.csv'.format(
            re.sub(r'\s+', '_', src_filename))
        self.ext = ext.lower()

        # Set directory as the same directory as the input file if not provided
        # For HDF5 files, create a further directory with the same name
        # as the input file since there will be multiple output files.
        self.dst = dst
        if self.dst is None:
            dst, _ = os.path.split(src)
            if self.ext.lower() == '.h5':
                dst = os.path.normpath(
                    os.path.join(os.path.expanduser(dst), self.src_filename))
            self.dst = dst
        else:
            self.dst = os.path.normpath(os.path.expanduser(dst))
        # Create directory tree if it does not exist and check for
        # read and write permissions.
        if not os.path.isdir(self.dst):
            logger.info("Creating directory '{}'".format(self.dst))
            os.makedirs(self.dst, exist_ok=True)
        logger.info("Checking write permission to directory '{}'".format(
            self.dst))
        os.access(self.dst, os.W_OK)

        self.skip_header_rows = skip_header_rows
        self.skip_footer_rows = skip_footer_rows
        self.sheet_name = sheet_name
        self.score_column = score_column
        self.input_type = input_type
        self.one_based = one_based

        # Initialize sequence information.
        self._wt_sequence = None
        self._offset = None
        self.codons = None
        self.protein_sequence = None
        self.offset = offset
        self.wt_sequence = wt_sequence

    @property
    def wt_sequence(self):
        return self._wt_sequence

    @wt_sequence.setter
    def wt_sequence(self, value):
        if isinstance(value, tuple):
            seq, offset = value
            self.offset = offset
        else:
            seq = str(value).upper()
        # Initialize sequence information.
        if not constants.dna_re.fullmatch(seq):
            raise ValueError("{} is not a valid DNA sequence.".format(seq))
        self._wt_sequence = seq[self.offset:].upper()
        self.protein_sequence = utilities.translate_dna(self.wt_sequence)
        self.codons = list(utilities.slicer(self.wt_sequence, 3))

    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, offset):
        offset = int(offset)
        if offset < 0:
            raise ValueError("Offset must be not be negative.")
        self._offset = offset

    @property
    def extension(self):
        return self.ext.lower()

    @property
    def input_is_h5(self):
        return self.ext.lower() == '.h5'

    @property
    def input_is_scores_based(self):
        return self.input_type == constants.score_type

    @property
    def output_directory(self):
        return os.path.normpath(os.path.expanduser(self.dst))

    @property
    def output_file(self):
        return os.path.normpath(
            os.path.join(self.output_directory, self.dst_filename,))

    def convert(self):
        """
        Runs `parse_input` and saves the Mavedb-compliant result to file.
        """
        logger.info("Processing file {}".format(self.src))
        mave_df = self.parse_input(self.load_input_file())
        logger.info('Writing to {}'.format(self.output_file))
        mave_df.to_csv(self.output_file, sep=',', index=None, na_rep=np.NaN)

    def load_input_file(self):
        raise NotImplementedError()

    def parse_input(self, df):
        raise NotImplementedError()

    def parse_row(self, row):
        raise NotImplementedError()

    def validate_against_wt_sequence(self, variant):
        """
        Checks that the reference base in a substitution variant matches that
        in the wild-type sequence provided.

        Parameters
        ----------
        variant : str
            A nucleotide substitution variant with valid HGVS_ syntax.
        """
        if is_multi(variant):
            _ = [
                self.validate_against_wt_sequence(v)
                for v in utilities.split_variant(variant)
            ]
            return

        if variant in constants.special_variants:
            return

        variant = utilities.NucleotideSubstitutionEvent(variant)
        if variant.silent:
            return

        zero_based_pos = variant.position - int(self.one_based)
        if zero_based_pos < 0:
            raise IndexError(
                "Encountered a negative position. "
                "Positions might not be one-based.")

        wt_ref_nt = self.wt_sequence[zero_based_pos]
        if variant.ref != wt_ref_nt:
            raise ValueError(
                "Base '{base}' at 1-based position {pos} in the wild-type "
                "sequence (offset {offset}) does not match the base "
                "suggested in the variant '{variant}'.".format(
                    pos=zero_based_pos + 1, base=wt_ref_nt, offset=self.offset,
                    variant=variant))

    def validate_against_protein_sequence(self, variant):
        """
        Checks that the reference amino acid in a substitution variant matches
        that in the translated wild-type sequence provided.

        Parameters
        ----------
        variant : str
            A protein substitution variant with valid HGVS_ syntax.
        """
        if is_multi(variant):
            _ = [
                self.validate_against_protein_sequence(v)
                for v in utilities.split_variant(variant)
            ]
            return

        if variant in constants.special_variants or 'p.=' in variant:
            return

        variant = utilities.ProteinSubstitutionEvent(variant)
        zero_based_pos = variant.position - int(self.one_based)
        if zero_based_pos < 0:
            raise IndexError(
                "Encountered a negative position. "
                "Positions might not be one-based.")

        wt_aa = constants.AA_CODES[self.protein_sequence[zero_based_pos]]
        if variant.ref != wt_aa:
            raise ValueError(
                "Amino acid '{aa}' at 1-based position {pos} in the "
                "translated protein sequence does not match the amino acid "
                "suggested in the variant '{variant}'.".format(
                    pos=zero_based_pos + 1, aa=wt_aa, variant=variant))
