import logging

import numpy as np

from . import LOGGER, constants, utilities


logger = logging.getLogger(LOGGER)


def drop_na_columns(df, inplace=False):
    """
    Drop columns where all entries are null. Operation is performed in place.
    """
    if not inplace:
        df = utilities.copy_dataframe(df)

    has_nt_col = constants.hgvs_nt_col in df.columns
    has_pro_col = constants.hgvs_pro_col in df.columns

    if has_nt_col:
        nt_all_null = np.all(
            df.loc[:, constants.hgvs_nt_col].apply(variant_utils.is_null))
        if nt_all_null:
            df.drop(columns=[constants.hgvs_nt_col], inplace=True)
    if has_pro_col:
        pro_all_null = np.all(
            df.loc[:, constants.hgvs_pro_col].apply(variant_utils.is_null))
        if pro_all_null:
            df.drop(columns=[constants.hgvs_pro_col], inplace=True)

    # Drop data columns that are all null.
    to_drop = list()
    for cname in utilities.non_hgvs_columns(df.columns):
        if np.all(df.loc[:, cname].isnull()):
            logger.warning("Dropping column '{}' because it contains all null "
                           "values".format(cname))
            to_drop.append(cname)
    if len(to_drop) > 0:
        df.drop(columns=to_drop, inplace=True)

    return df


def drop_na_rows(df, inplace=False):
    """
    Drop rows where all non-HGVS entries are null. Operation is performed in
    place.
    """
    if not inplace:
        df = utilities.copy_dataframe(df)

    null_rows = df.loc[:, utilities.non_hgvs_columns(
        df.columns)].isnull().all(axis=1)
    if sum(null_rows) > 0:
        logger.warning("Dropping {} rows that contain all null values".format(
            sum(null_rows)))
        df.drop(index=df.index[null_rows], inplace=True)

    return df
