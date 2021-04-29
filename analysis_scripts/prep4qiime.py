#!/usr/bin/env python3

import pandas as pd
import argparse
import os 

def extract_taxa(asv):
    taxa = asv[['taxonomy']]
    asv = asv.drop('taxonomy', axis=1)
    asv_out = asv.T
    return asv_out, taxa 

def main():
    """
    parse cli arguments and write output to file
    :return:
    """

    parser = argparse.ArgumentParser(description="DESCRIPTION\n"
                                                "Prepare DADA2 output for import into qiime\n"
                                                "\n\n==========================BASIC USAGE==========================\n"
                                                "\n$ prep4qiime.py --asv_table ASV_seqtab.tab -output qiime2 \n"
                                                , formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-a", "--asv_table", type=str, help="asv table path")
    parser.add_argument("-o", "--output", type=str, help="output folder")

    args = parser.parse_args()
    
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    asv_in = pd.read_csv(args.asv_table,sep="\t", header=0, index_col=0)

    asv_out, taxa_out = extract_taxa(asv_in)

    asv_outname = os.path.join(args.output, "qiime2_in.tsv")    
    taxa_outname = os.path.join(args.output, "qiime2_taxa.tsv") 
    asv_out.to_csv(asv_outname, sep="\t")
    taxa_out.to_csv(taxa_outname, sep="\t")


if __name__ == '__main__':
    main()
