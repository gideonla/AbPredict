#!/usr/bin/env python3


import pandas as pd 


def read_score_files(paths: list, index_col='description'):
    """
    :param paths: list of paths of score files
    :param index_col: the column to be used as a row name 
    :return: dataframe of the score files. 
    """
    frames = [pd.read_csv(sf, header=1, sep='\s+', index_col=index_col) 
              for sf in paths]
    return pd.concat(frames)
