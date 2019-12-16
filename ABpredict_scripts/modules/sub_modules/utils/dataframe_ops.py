#! /usr/bin/env python3
"""
commonly used list related operations
"""
import numpy as np
import pandas as pd


def sort_df(df_: pd.DataFrame, column_idx: str, key) -> pd.DataFrame:
    '''Takes dataframe, column index and custom function for sorting,
    returns dataframe sorted by this column using this function'''

    col = df_.ix[:, column_idx]
    temp = np.array(col.values.tolist())
    order = sorted(range(len(temp)), key=lambda j: key(temp[j]))
    return df_.ix[order]
