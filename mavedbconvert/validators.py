import logging

import hgvsp

from hgvs.sequencevariant import SequenceVariant
from hgvs.exceptions import HGVSParseError

import numpy as np
import pandas as pd
from numpy.testing import assert_array_equal

from tqdm import tqdm

from joblib import Parallel, delayed

from . import constants, utilities, exceptions, LOGGER


logger = logging.getLogger(LOGGER)


class ValidationBackend(object):
    """
    Validation backend which provides the interface `validate` for validating
    HGVS_ variants.
    """

    def validate(self, variant):
        raise NotImplementedError()


class HGVSPatternsBackend(ValidationBackend):
    """
    Backend using the regex based validation in `hgvsp`. Fast but may be
    too strict in specific cases.
    """

    def validate(self, variant):
        """
        Validates a HGVS_ variant using the regex patterns in `hgvsp`.

        Parameters
        ----------
        variant : str
            HGVS formatted variant string.

        Returns
        -------
        str
        """
        if variant in constants.special_variants:
            return variant
        single_match = hgvsp.single_variant_re.fullmatch(variant)
        multi_match = hgvsp.multi_variant_re.fullmatch(variant)
        if not (single_match or multi_match):
            raise exceptions.HGVSValidationError(
                "'{}' is not valid HGVS syntax.".format(variant))
        return variant


class HGVSBiocommonsBackend(ValidationBackend):
    """
    Backend using the grammar based validation in `hgvs`. Can be slow but is
    more robust.
    """

    def __init__(self, transcript=None):
        self.transcript = transcript or constants.dummy_ref

    def validate(self, variant):
        """
        Splits a variant if it is multi-variant to validate each individual
        variant since `hgvs` does not support multi syntax. Validates each
        HGVS variant against a set grammar.

        Parameters
        ----------
        variant : str
            HGVS formatted variant string.

        Returns
        -------
        list[`SequenceVariant`]
            List of sequence variants. Singular list if variant was not
            in multi-variant syntax.
        """
        if variant in constants.special_variants:
            return variant
        try:
            seqvars = []
            for v in utilities.split_variant(variant):
                seqvar = constants.hgvs_parser.parse_hgvs_variant(
                    '{}:{}'.format(self.transcript, v))
                seqvar.validate()
                seqvars.append(seqvar)
            return seqvars if len(seqvars) > 1 else seqvars[0]
        except HGVSParseError as e:
            raise exceptions.HGVSValidationError(
                "'{}' is not valid HGVS syntax "
                "for the following reason: {}".format(variant, e))


def validate_variants(variants, transcript=None, validation_backend=None,
                      n_jobs=1, verbose=0, backend='multiprocessing', ):
    """
    Validate each variant's HGVS_ syntax.

    Parameters
    ----------
    variants : list[str]
        Variant HGVS_ representations.
    transcript : str, optional.
        Transcript the variants reference.
    validation_backend : ValidationBackend
        A parsing backend implementing `validate`.
    n_jobs : int, optional
        Number of jobs to run in parallel.
    verbose : int, optional
        Joblib's verbosity level.
    backend : str, optional
        Parallel backend to use. Defaults to `multiprocessing`.

    Returns
    -------
    list[Union[str, SequenceVariant]]
        Formatted and validated variants.
    """
    if validation_backend is None:
        if transcript is None:
            validation_backend = HGVSPatternsBackend()
        else:
            validation_backend = HGVSBiocommonsBackend(transcript)
    return Parallel(n_jobs=n_jobs, verbose=verbose, backend=backend)(
        delayed(validation_backend.validate)(variant)
        for variant in variants
    )


def validate_has_column(df, column):
    """Validates that a `DataFrame` contains `column` in it's columns."""
    if column not in df.columns:
        raise KeyError(
            "Missing column '{}'. Existing columns are {}.".format(
                column, ', '.join(df.columns)
            )
        )


def validate_columns_are_numeric(df):
    """Checks non-hgvs columns for float or int data."""
    for column in df.columns:
        if column in [constants.hgvs_pro_col, constants.hgvs_nt_col]:
            continue
        else:
            if not (np.issubdtype(df.dtypes[column], np.floating) or
                    np.issubdtype(df.dtypes[column], np.integer)):
                raise TypeError(
                    "Expected only float or int data columns. Got {}.".format(
                        str(df.dtypes[column])
                    ))


def validate_hgvs_nt_uniqueness(df):
    """Validate that hgvs columns only define a variant once."""
    df = df[~df[constants.hgvs_nt_col].isnull()]
    dups = df.loc[:, constants.hgvs_nt_col].duplicated(keep=False)
    if np.any(dups):
        dup_list = ["{} ({})".format(x, y) for x, y in
                    zip(df.loc[dups, constants.hgvs_nt_col],
                        dups.index[dups])]
        raise ValueError(
            "duplicate HGVS nucleotide strings found: {}".format(
                ', '.join(sorted(dup_list))
            ))


def validate_hgvs_pro_uniqueness(df):
    """Validate that hgvs columns only define a variant once."""
    df = df[~df[constants.hgvs_pro_col].isnull()]
    dups = df.loc[:, constants.hgvs_pro_col].duplicated(keep=False)
    if np.any(dups):
        dup_list = ["{} ({})".format(x, y) for x, y in
                    zip(df.loc[dups, constants.hgvs_pro_col],
                        dups.index[dups])]
        raise ValueError(
            "Duplicate HGVS protein strings found: {}".format(
                ', '.join(sorted(dup_list))
            ))


def validate_datasets_define_same_variants(scores_df, counts_df):
    """
    Checks if two `pd.DataFrame` objects parsed from uploaded files
    define the same variants.

    Parameters
    ----------
    scores_df : `pd.DataFrame`
        Scores dataframe parsed from an uploaded scores file.
    counts_df : `pd.DataFrame`
        Scores dataframe parsed from an uploaded counts file.
    """

    scores_columns = [c for c in scores_df.columns
                      if c in constants.variant_columns]
    counts_columns = [c for c in counts_df.columns
                      if c in constants.variant_columns]
    if scores_columns != counts_columns:
        raise AssertionError(
            "Dataframes define different hgvs columns. "
            "Scores defines '{}' and counts defines '{}'.".format(
                ', '.join(scores_columns), ', '.join(counts_columns)
            )
        )

    if constants.hgvs_nt_col in scores_columns:
        scores_nt = scores_df[constants.hgvs_nt_col].values
        counts_nt = counts_df[constants.hgvs_nt_col].values
        try:
            assert_array_equal(scores_nt, counts_nt)  # Treats np.NaN as equal
        except AssertionError:
            not_equal_selector = scores_nt != counts_nt
            neq_list = [
                "{} ({})".format(x, y) for x, y in
                zip(scores_nt[not_equal_selector],
                    counts_nt[not_equal_selector])
                if (x is not np.NaN) and (y is not np.NaN)
            ]
            raise AssertionError(
                "Scores and counts do not define the same "
                "nucleotide variants: {}.".format(', '.join(neq_list))
            )

    if constants.hgvs_pro_col in scores_columns:
        scores_pro = scores_df[constants.hgvs_pro_col].values
        counts_pro = counts_df[constants.hgvs_pro_col].values
        try:
            assert_array_equal(scores_pro,
                               counts_pro)  # Treats np.NaN as equal
        except AssertionError:
            not_equal_selector = scores_pro != counts_pro
            neq_list = [
                "{} ({})".format(x, y) for x, y in
                zip(scores_pro[not_equal_selector],
                    counts_pro[not_equal_selector])
                if (x is not np.NaN) and (y is not np.NaN)
            ]
            raise AssertionError(
                "Scores and counts do not define the same protein variants: "
                "{}.".format(', '.join(neq_list))
            )


def validate_mavedb_compliance(df, df_type):
    """Runs MaveDB compliance checks."""
    tqdm.pandas(desc="Validating variants")

    has_nt_col = constants.hgvs_nt_col in df.columns
    has_pro_col = constants.hgvs_pro_col in df.columns
    if not has_nt_col and not has_pro_col:
        raise ValueError(
            "Dataframe must define either '{}', '{}' or both.".format(
                constants.hgvs_nt_col, constants.hgvs_pro_col
            ))

    primary_col = None
    if has_nt_col:
        defines_nt = not all(
            df.loc[:, constants.hgvs_nt_col].progress_apply(utilities.is_null))
        if defines_nt:
            primary_col = constants.hgvs_nt_col

    if has_pro_col and primary_col is None:
        defines_pro = not all(
            df.loc[:, constants.hgvs_pro_col].progress_apply(
                utilities.is_null))
        if defines_pro:
            primary_col = constants.hgvs_pro_col

    if primary_col is None:
        raise ValueError("Neither '{}' or '{}' defined any variants.".format(
            constants.hgvs_nt_col, constants.hgvs_pro_col
        ))

    null_primary = df.loc[:, primary_col].progress_apply(utilities.is_null)
    if any(null_primary):
        raise ValueError(
            "Primary column (inferred as '{}') cannot "
            "contain the null values {} (case-insensitive).".format(
                primary_col, 'NaN, Na, None, whitespace, Undefined'
            ))

    if primary_col == constants.hgvs_nt_col:
        validate_hgvs_nt_uniqueness(df)
    elif primary_col == constants.hgvs_pro_col:
        try:
            validate_hgvs_pro_uniqueness(df)
        except ValueError as e:
            logger.warning(
                "'{}' column contains duplicated entries: {}".format(
                    constants.hgvs_pro_col, e))

    validate_columns_are_numeric(df)
    if df_type == constants.score_type:
        validate_has_column(df, 'score')
    return df
