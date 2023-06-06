#!/usr/bin/env python

import argparse
import pandas as pd


## main function
# script needs to load in sequence stats files pre and post-QC and combine to generate final report 
def summarise_by_sample(df):
    # add sample name
    df['sample'] = df['file'].str.split('_').str[0]
    df_summarised = df.groupby(['sample'])[['num_seqs','sum_len']].sum()
    # divide sum_len by 1e9
    df_summarised['sum_len'] = round(df_summarised['sum_len'] / 1e9, 2)
    return df_summarised

def create_report(raw, trimmed, clean):

    raw_stats_sum  =  summarise_by_sample(raw)
    trimming_stats_sum = summarise_by_sample(trimmed)
    clean_stats_sum = summarise_by_sample(clean)
    
    raw_stats_sum.rename(columns={'num_seqs': 'num_seqs_raw', 'sum_len': 'num_bases_raw'}, inplace=True)
    trimming_stats_sum.rename(columns={'num_seqs': 'num_seqs_trimmed', 'sum_len': 'num_bases_trimmed'}, inplace=True)
    clean_stats_sum.rename(columns={'num_seqs': 'num_seqs_clean', 'sum_len': 'num_bases_clean'}, inplace=True)

    out = pd.concat([raw_stats_sum, trimming_stats_sum, clean_stats_sum], axis=1)
    return out 

def main():
    """
    parse cli arguments and write output to file
    :return:
    """
    parser = argparse.ArgumentParser(description="DESCRIPTION\n"
                                    "This script takes in multiple seqkit stats reports for the various steps of shotgun\n" 
                                    "sequences and generates a combined report.\n"
                                    "\n\n==========================BASIC USAGE==========================\n"
                                    "\n$ sequence_stats.py -r sequence_stats_raw.tsv -t sequence_stats_trimmed.tsv\n"
                                    "-c sequence_stats_clean.tsv -o sequence_stats.tsv", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-r", "--raw", required=True, type=str, help="seqkit report on raw reads")
    parser.add_argument("-t", "--trimmed", required=True, type=str, help="seqkit report on trimmed reads")
    parser.add_argument("-c", "--cleaned", required=True, type=str, help="seqkit report on cleaned reads")
    parser.add_argument("-o", "--output", required=True, type=str, help="output file name")
    args = parser.parse_args()
    
    raw_stats = pd.read_csv(args.raw, sep="\t")
    trimming_stats = pd.read_csv(args.trimmed, sep="\t")
    clean_stats = pd.read_csv(args.cleaned, sep="\t")

    combined_report = create_report(raw_stats, trimming_stats, clean_stats)

    combined_report.to_csv(args.output, sep="\t")

if __name__ == '__main__':
    main()