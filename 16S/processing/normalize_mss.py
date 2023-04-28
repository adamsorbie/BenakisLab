#!/usr/bin/env python3

import qiime2
from qiime2.plugins import feature_table
import pandas as pd
import sys
import numpy as np


# why does this transpose? 
# get identifier from unzipped qiime2 artifact
def normalize(otu, rel=False):
    """
    :param  pandas.Dataframe otu: OTU table to be normalised
    :param rel: Whether or not to return OTU normalised by relative abundance
    :return: OTU table normalised by minimum sum scaling or relative abundance
    """
    if otu.shape[1] == otu.select_dtypes(include=np.number).shape[1]:
        sum_row = otu.sum()
        if rel:
            rel_abund = otu * 100 / sum_row
            return rel_abund
        else:
            min_sum = min(sum_row)
            normalised = otu / sum_row * min_sum
            return normalised
    else:
        print("Non-numeric values are not allowed")


data_qza = qiime2.Artifact.load(sys.argv[1])
# transpose so after loading artifact and converting to df, samples are columns and features rows 
df = data_qza.view(pd.DataFrame).T
normalized_df = normalize(df)
# get colSums
colsums = normalized_df.sum(axis=0)


normalized_df.to_csv("feature-table-norm.tsv", sep="\t")
colsums.to_csv('logfile.tsv', sep="\t")

