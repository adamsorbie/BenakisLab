#!/usr/bin/env python3
import pandas as pd
import numpy as np 
import sys 

# Functions 

## normalize data by converting to relative abundance
def normalize(asvtab):
    sum_row = asvtab.sum()
    rel_abund = asvtab * 100 / sum_row
    return rel_abund
# filter by relative abundance
def filter_abund(asvtab, abun_thresh=0.25):
    asvtab_filt = asvtab[asvtab.sum(axis=1) > abun_thresh]
    asv_idx = asvtab_filt.index.tolist()
    return(asv_idx)
# filter asv table and return filtered table and list of ASVs to keep 
def threshold_asvtab(asvtab):
    asvtab_rel = normalize(asvtab)
    filt_asvs = filter_abund(asvtab_rel)
    asvtab_out = asvtab[asvtab.index.isin(filt_asvs)]
    return(asvtab_out, filt_asvs)

# read in file (with annoying biom header removed)
asvtable = pd.read_csv("feature-table-2.tsv", sep="\t", index_col=0)
# filter asv table 
asvtab_filt, asv_filt_list = threshold_asvtab(asvtable)
# write out filtered file  
asvtab_filt.to_csv("feature-table-filt.tsv", sep="\t")
# write out list of ASVs for filtering fasta
asv_filt_list_out = pd.DataFrame(asv_filt_list)
asv_filt_list_out[0] = '>' + asv_filt_list_out[0].astype(str)
asv_filt_list_out.to_csv("keep-asvs.txt", index=False, sep="\t", header=False)




