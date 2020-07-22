import logging
from tqdm import tqdm

import pandas as pd
import numpy as np
from fqfa.constants.iupac.protein import AA_CODES


from . import LOGGER, constants, base, utilities, filters, validators


__all__ = ["Enrich"]


logger = logging.getLogger(LOGGER)


class Enrich(base.BaseProgram):
    __doc__ = base.BaseProgram.__doc__

    def __init__(
        self,
        src,
        wt_sequence,
        offset=0,
        dst=None,
        one_based=False,
        skip_header_rows=0,
        skip_footer_rows=0,
        score_column=None,
        hgvs_column=None,
        input_type=None,
        sheet_name=None,
        is_coding=True,
    ):
        super().__init__(
            src=src,
            wt_sequence=wt_sequence,
            offset=offset,
            dst=dst,
            one_based=one_based,
            skip_header_rows=skip_header_rows,
            skip_footer_rows=skip_footer_rows,
            sheet_name=sheet_name,
            score_column=score_column,
            input_type=input_type,
            hgvs_column=hgvs_column,
            is_coding=is_coding,
        )
        if not abs(offset) % 3 == 0:
            raise ValueError("Enrich offset must be a multiple of 3.")

        if not is_coding:
            raise ValueError("Enrich does not support non-coding datasets.")

        if not self.score_column and self.input_type == constants.score_type:
            raise ValueError(
                "A score column must be specified if "
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
            logger.info("Skipping first {} row(s).".format(self.skip_footer_rows + 1))
        if self.skip_footer_rows:
            logger.info("Skipping last {} row(s).".format(self.skip_footer_rows + 1))

        if self.extension in (".xlsx", ".xls"):
            od = pd.read_excel(
                self.src,
                na_values=constants.extra_na,
                sheet_name=self.sheet_name,
                skiprows=self.skip_header_rows,
                skipfooter=self.skip_footer_rows,
            )
            if not self.sheet_name:
                self.sheet_name = list(od.keys())[0]
                if len(od) > 1:
                    logger.warning(
                        "Multiple sheet names detected ({}). Parsing "
                        "{} only. Re-run with the `sheet_name` argument "
                        "to parse a specific sheet.".format(
                            ", ".join(list(od.keys())), self.sheet_name
                        )
                    )
            df = od[self.sheet_name]
        else:
            sep = "\t"
            if self.ext.lower() == ".csv":
                sep = ","
            df = pd.read_csv(
                self.src,
                delimiter=sep,
                na_values=constants.extra_na,
                skipfooter=self.skip_footer_rows,
                skiprows=self.skip_header_rows,
            )

        if "seqID" not in df.columns:
            raise ValueError("Input is missing the required column 'seqID'.")
        return df

    def parse_row(self, row):
        """
        Parse an Enrich seq_id value in the format:
            `<comma delim position list>-<comma delim aa list>`

        Into valid HGVS formated nucleotide and protein strings.

        Parameters
        ----------
        row :  str
            `Enrich` formatted SeqID value.

        Returns
        -------
        `str`
            An hgvs_pro string. The positions reported will be 1-based.
        """
        seq_id = row
        if not seq_id or "NA" in seq_id.upper():
            raise ValueError("'{}' is a malformed SeqID.".format(seq_id))

        positions, aa_codes = seq_id.split("-")
        positions = positions.split(",")
        aa_codes = aa_codes.split(",")
        events = []

        if len(positions) != len(aa_codes):
            raise ValueError(
                "Number of positions ({pos}) in {seqid} "
                "does not match number of ammino acid codes ({codes}).".format(
                    pos=len(positions), seqid=seq_id, codes=len(aa_codes)
                )
            )

        for position, aa in zip(positions, aa_codes):
            offset = (1, -1)[self.offset < 0] * abs(self.offset) // 3
            aa_position = int(position) - int(self.one_based) + 1 - offset
            if aa_position < 1:
                raise IndexError(
                    "Position in SeqID '{pos}-{aa}' from row '{seqid}' must "
                    "be 1 or greater after applying codon adjusted offset "
                    "{offset} ({raw_offset} / 3). Computed position is "
                    "{aa_pos}.".format(
                        pos=position,
                        aa=aa,
                        seqid=seq_id,
                        raw_offset=self.offset,
                        offset=offset,
                        aa_pos=aa_position,
                    )
                )
            if aa_position > len(self.protein_sequence):
                raise IndexError(
                    "Position in SeqID '{pos}-{aa}' from row '{seqid}' is "
                    "out of bounds after applying codon adjusted offset "
                    "{offset} ({raw_offset} / 3). Computed position is "
                    "{aa_pos} and the length of the translated sequence "
                    "is {seqlen}.".format(
                        pos=position,
                        aa=aa,
                        seqid=seq_id,
                        raw_offset=self.offset,
                        offset=offset,
                        aa_pos=aa_position,
                        seqlen=len(self.protein_sequence),
                    )
                )

            wt_aa = AA_CODES[self.protein_sequence[aa_position - 1].upper()]
            if aa == "?":
                mut_aa = "???"
            else:
                mut_aa = AA_CODES[aa.upper()]
            if wt_aa == mut_aa:
                events.append("{wt}{pos}=".format(wt=wt_aa, pos=aa_position))
            else:
                events.append(
                    "{wt}{pos}{mut}".format(wt=wt_aa, pos=aa_position, mut=mut_aa)
                )

        if len(events) == 0:
            raise ValueError(
                "Could not parse any variant strings from {}".format(seq_id)
            )
        return utilities.hgvs_pro_from_event_list(events)

    def parse_input(self, df):
        """
        Parse a list of Enrich seq_id values in the format:
            `<comma delim position list>-<comma delim aa list>`

        Into a list of valid HGVS formated nucleotide and protein strings.

        Parameters
        ----------
        df :  pd.DataFrame
            DataFrame containing data in Enrich format. Expects that the
            DataFrame contains a column 'seqID' containing the seqIDs that will
            be converted to HGVS. All other columns are considered data columns
            and will be retained in the order that they appear, with the same
            column names.

        Returns
        -------
        `pd.DataFrame`
            A DataFrame containing the new HGVS strings and the previous data
            values suitable for MaveDB import.
        """
        data_columns = [c for c in df.columns if c != "seqID"]

        # output the conversion progress with a progress bar
        tqdm.pandas(desc="Parsing seqIDs")
        df.loc[:, constants.pro_variant_col] = df.loc[:, "seqID"].progress_apply(
            self.parse_row
        )

        # enrich output has no nucleotide data
        df.loc[:, constants.nt_variant_col] = None

        # reorder columns and drop seqID
        columns = list(constants.variant_columns)
        columns.extend(data_columns)
        # Duplicate of columns to remove entries during iteration if column
        # is non-numeric
        mave_columns = list(columns)
        data = {}
        for column in columns:
            column_type = df.dtypes[column]
            column_values = df[column].values

            if column in constants.variant_columns:
                astype = str
            elif np.issubdtype(column_type, np.floating):
                astype = np.float
            elif np.issubdtype(column_type, np.signedinteger):
                astype = np.int
            else:
                logger.warning("Dropping non-numeric column '{}'".format(column))
                mave_columns.remove(column)
                continue

            if self.input_is_scores_based and column == self.score_column:
                mave_columns.remove(column)
                column = constants.mavedb_score_column
            data[column] = list(utilities.format_column(column_values, astype))

        # Sort column order so 'score' comes right after hgvs columns.
        if self.input_is_scores_based:
            mave_columns = (
                mave_columns[:2] + [constants.mavedb_score_column] + mave_columns[2:]
            )
        mavedb_df = pd.DataFrame(data=data, columns=mave_columns)
        filters.drop_na_rows(mavedb_df)
        filters.drop_na_columns(mavedb_df)

        logger.info("Running MaveDB compliance validation.")
        validators.validate_mavedb_compliance(mavedb_df, df_type=self.input_type)
        return mavedb_df
