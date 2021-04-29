#!/usr/bin/env python3

import pandas as pd 
import sys 

otu_contam = pd.read_csv(sys.argv[1], header=0, index_col=0, sep="\t") 
noncontam_otus = pd.read_csv(sys.argv[2], header=None, index_col=0, sep="\t")

list_noncontam_otus = noncontam_otus.index.tolist()

otu_noncontam = otu_contam[otu_contam.index.isin(list_noncontam_otus)]

otu_noncontam.to_csv(sys.argv[3], sep="\t")

