#!/usr/bin/env python

import pandas as pd
import sys
import numpy as np

def cat_id_description(df, filetype):
    if filetype == "MG":
        df = df.reset_index()
        df['function-description'] = df['function'].map(str) + "-" + df['description']
        df = df.drop(['function', 'description'], axis=1)
        df.set_index('function-description', inplace=True)
        return df
    elif filetype == 'PWY':
        df = df.reset_index()
        df['pathway-description'] = df['pathway'].map(str) + "-" + df['description']
        df = df.drop(['pathway', 'description'], axis=1)
        df.set_index('pathway-description', inplace=True)
        return df
    else:
        print("ERROR, NO FILETYPE SPECIFIED")

# add taxonomy column
def add_taxonomy(pred_meta, taxonomy):
    # make taxonomy dict mapping otus to taxonomy string
    taxmap =  pd.Series(taxonomy.taxonomy.values,index=taxonomy.index).to_dict()
    pred_meta['taxonomy'] = pred_meta['taxon'].map(taxmap)
    return pred_meta

def format_taxonomy(taxonomy):
    # split into columns for each rank
    taxonomy[["Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"]] = taxonomy['taxonomy'].str.split(';',expand=True)
	# replace empty cells with nan
    taxonomy.replace("", np.nan, inplace=True)
	# get empty indices
    taxonomy_empty = taxonomy.isnull()
	# forward fill na
    taxonomy_format =  taxonomy.fillna(method='ffill', axis=1)
	# using indices from before add prefix to ffilled cells
    taxonomy_format.update('unknown_' + taxonomy_format.mask(~taxonomy_empty))
	# combine formatted taxonomy into one new column, overwriting old taxonomy
    taxonomy_format["taxonomy"] = taxonomy_format[taxonomy_format.columns[1:]].apply(lambda x: ";".join(x.dropna().astype(str)),
    axis=1)
	# drop extra cols and return
    taxonomy_format.drop(["Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"], axis=1, inplace=True)

    return taxonomy_format

# format all at once
def format_for_r(pred_meta, taxonomy, filetype, strat="unstrat"):
    if strat == "strat":
        pred_meta_in = add_taxonomy(pred_meta, taxonomy)
    elif strat == "unstrat":
        pred_meta_in = pred_meta
    else:
        print("Error please enter valid value. \n Options:"
              "strat, unstrat")
    pred_meta_out = cat_id_description(pred_meta_in, filetype=filetype)
    return pred_meta_out


# read files
pred_meta = pd.read_csv(sys.argv[1], sep="\t", header=0, index_col=0)
taxonomy = pd.read_csv(sys.argv[2], sep="\t",
                      header=0, index_col=0, usecols=['Unnamed: 0', "taxonomy"])

taxonomy_format = format_taxonomy(taxonomy)

pred_meta_out = format_for_r(pred_meta, taxonomy_format, filetype=sys.argv[3], strat=sys.argv[4])

pred_meta_out.to_csv("pred_metagenome_taxa.tsv", sep="\t")
