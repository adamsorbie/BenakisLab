# this script removes negatives from an OTU table using a sample mapping file 

import pandas as pd 
import argparse

# read in otu table and mapping file
def read_files(metadata, otu):
        metadata = pd.read_csv(metadata, index_col=0, header=0, sep="\t")
        otu = pd.read_csv(otu, index_col=0, header=0, sep="\t")
        return metadata, otu

def main():
        """
        parse cli arguments and write output to file
        :return:
        """
        parser = argparse.ArgumentParser(description="DESCRIPTION\n"
                                        "This script automatically removes negative controls \n"
                                        "and low read samples based on a user provided threshold and  \n"
                                        "mapping file containing the samples for analysis\n"
                                        "\n\n==========================BASIC USAGE==========================\n"
                                        "\n$ remove_negatives.py -i otu.tab -m meta.tab -o otu_neg_rem.tab\n"
                                        ,formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("-i", "--input", required=True, type=str, help="OTU file path")
        parser.add_argument("-m", "--mapping", required=True, type=str, help="Mapping file path")
        parser.add_argument("-o", "--output", required=False, type=str, help="output file name")
        parser.add_argument("-r", "--read_threshold", required=False, type=int, help="output file name")
        args = parser.parse_args()
        
        meta, otu_tab = read_files(args.mapping, args.input)
        
        tax_col = otu_tab[["taxonomy"]]
        
        samples = meta.index.tolist()
        
        otu_neg_rem = otu_tab[samples]
        
        if args.read_threshold:
                otu_low_reads_rem = otu_neg_rem[otu_neg_rem.columns[otu_neg_rem.sum() > args.read_threshold]]
                keep_samples = otu_low_reads_rem.columns.tolist()
                # also remove low read samples from mapping file
                meta_filt = meta[meta.index.isin(keep_samples)]
                meta_filt.to_csv("mapping_low_reads_rem.tab", sep="\t")
        else:
                otu_low_reads_rem = otu_neg_rem
        
        otu_low_reads_rem["taxonomy"] = tax_col
        
        if args.output == None:
                otu_low_reads_rem.to_csv("OTUs-Table_wo_ctrls.tab", sep="\t")
        else:
                otu_low_reads_rem.to_csv(args.output, sep="\t") 
        


if __name__ == '__main__':
        main()