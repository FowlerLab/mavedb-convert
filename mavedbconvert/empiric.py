import logging
from tqdm import tqdm

import pandas as pd
import numpy as np

from . import base, utilities, constants, filters, validators, LOGGER


logger = logging.getLogger(LOGGER)


__all__ = [
    'Empiric',
    'infer_nt_substitution',
    'infer_pro_substitution',
]


def infer_nt_substitution(wt_codon, mut_codon, codon_pos):
    """
    Following block will report substitutiom events by comparing
    nucleotides within the wt and mutant codons. The position of the
    nucleotide level event is computed as:

        3 * 0-based codon position + 1-based pos of nt within codon

    Parameters
    ----------
    wt_codon : `str`
        The codon from the wild-type sequence.

    mut_codon : `str`
        The codon from the mutant-type sequence.

    codon_pos : `int`
        The 0-based position of the codon in the wild-type sequence.

    Returns
    -------
    `str`
        The inferred coding DNA HGVS-formatted string.
    """
    events = []
    for i, nt in enumerate(wt_codon):
        pos = (3 * codon_pos) + (i + 1)
        if nt.lower() != mut_codon[i].lower():
            events.append(
                '{pos}{wt_nt}>{mut_nt}'.format(
                    wt_nt=nt.upper(), pos=pos, mut_nt=mut_codon[i].upper())
            )
        else:
            events.append('{pos}='.format(pos=pos))
    return utilities.hgvs_nt_from_event_list(events, prefix='c')


def infer_pro_substitution(wt_aa, mut_aa, codon_pos):
    """
    Infer a HGVS-formatted subsitution event based on the wild-type amino acid
    and mutant-type amino acid.

    Parameters
    ----------
    wt_aa : `str`
        The amino acid from the wild-type protein sequence.

    mut_aa : `str`
        The amino acid from the mutant-type protein sequence.

    codon_pos : `int`
        The 0-based position of the codon in the wild-type sequence. This
        is equal to the position of the amino acid in the wild-type sequence.

    Returns
    -------
    `str`
        The HGVS-formatted subsitution event.
    """
    wt_aa = constants.AA_CODES[wt_aa.upper()]
    mut_aa = constants.AA_CODES[mut_aa.upper()]
    if wt_aa.lower() == mut_aa.lower():
        return utilities.hgvs_pro_from_event_list(['{wt_aa}{pos}='.format(
            wt_aa=wt_aa, pos=codon_pos + 1)])
    else:
        return utilities.hgvs_pro_from_event_list(
            ['{wt_aa}{pos}{mut_aa}'.format(
                wt_aa=wt_aa, pos=codon_pos + 1, mut_aa=mut_aa)]
        )


class Empiric(base.BaseProgram):
    __doc__ = base.BaseProgram.__doc__

    def __init__(self, src, wt_sequence, dst=None, one_based=False,
                 skip_header_rows=0, skip_footer_rows=0, score_column=None,
                 input_type=None, sheet_name=None, **kwargs):
        super().__init__(
            src=src,
            wt_sequence=wt_sequence,
            offset=0,
            dst=dst,
            one_based=one_based,
            skip_header_rows=skip_header_rows,
            skip_footer_rows=skip_footer_rows,
            sheet_name=sheet_name,
            score_column=score_column,
            input_type=input_type,
        )
        if not self.score_column and self.input_type == constants.score_type:
            raise ValueError(
                "A score column must be sepcified if "
                "the input file is a scores file."
            )

    def load_input_file(self):
        """
        Loads the input file specified at initialization into a dataframe.

        Returns
        -------
        `pd.DataFrame`
        """
        if self.skip_header_rows:
            logger.info("Skipping first {} row(s).".format(
                self.skip_footer_rows + 1))
        if self.skip_footer_rows:
            logger.info("Skipping last {} row(s).".format(
                self.skip_footer_rows + 1))
            
        if self.extension in ('.xlsx', '.xls'):
            od = pd.read_excel(
                self.src,
                na_values=constants.extra_na,
                skiprows=self.skip_header_rows,
                skipfooter=self.skip_footer_rows,
                sheet_name=self.sheet_name
            )
            if not self.sheet_name:
                self.sheet_name = list(od.keys())[0]
                if len(od) > 1:
                    logger.warning(
                        "Multiple sheet names detected ({}). Parsing "
                        "{} only. Re-run with the `sheet_name` argument "
                        "to parse a specific sheet.".format(
                            ', '.join(list(od.keys())), self.sheet_name)
                    )
            return od[self.sheet_name]
        else:
            df = pd.read_csv(
                self.src,
                delimiter='\t',
                na_values=constants.extra_na,
                skipfooter=self.skip_footer_rows,
                skiprows=self.skip_header_rows,
            )
        if 'Position' not in df.columns:
            raise ValueError(
                "Input is missing the required column 'Position'.")
        if 'Amino Acid' not in df.columns:
            raise ValueError(
                "Input is missing the required column 'Amino Acid'.")
        if 'Codon' not in df.columns:
            logger.warning(
                "Warning: Input is missing the column 'Codon'. "
                "Nucelotide level variants will not be inferred."
            )
        return df

    def parse_row(self, row):
        """
        Parses a dataframe row containing the columns 'Position', 'Amino Acid'
        and optionally 'Codon'.

        Parameters
        ----------
        row : `pd.DataFrame`
            A row from a pd.DataFrame

        Returns
        -------
        `tuple`
            A 2-tuple (hgvs_nt, hgvs_pro), where hgvs_nt will be `None` if
            `infer_nt` is `False`.
        """
        infer_nt = 'Codon' in row.axes[0]
        codon_pos = int(row['Position']) - int(self.one_based)
        mut_aa = str(row['Amino Acid']).strip().upper()
        if codon_pos < 0:
            raise IndexError(
                "Negative position encountered after adjusting for 1-based "
                "coordinate input. Coordinates might not be one-based."
            )
        if codon_pos > len(self.codons) - 1:
            raise IndexError(
                "Coordinate {pos} (1-based) is out of bounds. The maximum "
                "index of the translated sequence is {idx} "
                "(length {length}).".format(
                    pos=codon_pos + 1,
                    idx=len(self.codons),
                    length=len(self.codons)))

        if utilities.is_null(mut_aa) or not mut_aa:
            raise ValueError(
                "Missing `Amino Acid` value in row '{}'.".format(
                    row['row_num'])
            )

        wt_codon = self.codons[codon_pos].upper()
        wt_aa = constants.CODON_TABLE[wt_codon].upper()
        if infer_nt:
            mut_codon = str(row['Codon']).strip().upper()
            if utilities.is_null(mut_codon) or not mut_codon:
                raise ValueError(
                    "Missing `Codon` value in row '{}'.".format(row['row_num'])
                )
            if constants.CODON_TABLE[mut_codon] != mut_aa:
                raise ValueError(
                    "Codon '{}' does not match the amino acid '{}' "
                    "specified in the 'Amino Acid' column in row {}.".format(
                        mut_codon, mut_aa, row['row_num']))
            hgvs_nt = infer_nt_substitution(wt_codon, mut_codon, codon_pos)
        else:
            hgvs_nt = None

        hgvs_pro = infer_pro_substitution(wt_aa, mut_aa, codon_pos)
        return hgvs_nt, hgvs_pro

    def parse_input(self, df):
        """
        Formats an input `pd.DataFrame` loaded from an `Empiric` formatted file
        into a format suitable for upload to `MaveDB`.

        Returns
        -------
        `pd.DataFrame`
            The `MaveDB` formatted dataframe.
        """
        df['row_num'] = range(0, len(df))

        rows = tqdm(
            df.iterrows(), desc='Parsing variants', total=len(df))
        tups = [self.parse_row(row) for _, row in rows]

        df[constants.nt_variant_col] = [tup[0] for tup in tups]
        df[constants.pro_variant_col] = [tup[1] for tup in tups]
        df.drop(columns=['Position', 'Amino Acid', 'row_num'], inplace=True)
        if 'Codon' in df.columns:
            df.drop(columns=['Codon'], inplace=True)

        data = {}
        mave_columns = list(constants.variant_columns) + \
            [c for c in df.columns if c not in constants.variant_columns]

        for column in df.columns:
            column_type = df.dtypes[column]
            column_values = df[column].values

            if column in constants.variant_columns:
                astype = str
            elif np.issubdtype(column_type, np.floating):
                astype = np.float
            elif np.issubdtype(column_type, np.signedinteger):
                astype = np.int
            else:
                logger.warning(
                    "Dropping non-numeric column '{}'".format(column))
                mave_columns.remove(column)
                continue

            if self.input_is_scores_based and column == self.score_column:
                mave_columns.remove(column)
                column = constants.mavedb_score_column
            data[column] = list(utilities.format_column(column_values, astype))

        # Sort column order so 'score' comes right after hgvs columns.
        if self.input_is_scores_based:
            mave_columns = mave_columns[:2] + \
                           [constants.mavedb_score_column] + \
                           mave_columns[2:]
        mavedb_df = pd.DataFrame(data=data, columns=mave_columns)
        filters.drop_na_rows(mavedb_df, inplace=True)
        filters.drop_na_columns(mavedb_df, inplace=True)

        logger.info("Running MaveDB compliance validation.")
        validators.validate_mavedb_compliance(
            mavedb_df, df_type=self.input_type)
        return mavedb_df
