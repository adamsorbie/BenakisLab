#!/bin/bash 
# Author: Adam Sorbie
# Date: 29/04/21
# Version: 0.5.2
 while getopts p:o: flag
 do
    case "${flag}" in
      p) path=${OPTARG};;
      o) out=${OPTARG};;
      *) echo "usage: $0 [-p] [-o]" >&2
        exit 1 ;;
    esac
 done

 if ((OPTIND == 1))
 then
    echo "No options specified, exiting"
    exit 1
 fi
 
eval "$(conda shell.bash hook)"
conda activate qiime2-2021.4
# strip trailing slash if existing 
outdir=$(echo ${out%/}) 

# declare variables - these can be hardcoded since pipeline output always has the same names
asv_tab=${path}"/ASV_seqtab_tax.tab"
asv_seqs=${path}"/ASV_seqs.fasta"
asv_seqs_aln=${path}"/aligned.fasta"
asv_tree=${path}/"ASV_tree.tre"

# grab taxonomy and split from otu table 
python prep4qiime.py -a $asv_tab -o $outdir

# cd to outdir since remaining commands should be in path 
cd $outdir || 
echo $PWD
biom convert -i qiime2_in.tsv -o qiime2_in.biom --table-type="OTU table" --to-hdf5 


## import files into qiime 

# import taxonomy 
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path qiime2_taxa.tsv \
--output-path taxonomy.qza

# import asv table 
qiime tools import \
--input-path qiime2_in.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path feature-table.qza

# import rep seqs 
qiime tools import \
--input-path $asv_seqs \
--type 'FeatureData[Sequence]' \
--output-path rep-seqs.qza

# import aligned seqs 
qiime tools import \
--input-path $asv_seqs_aln \
--type 'FeatureData[AlignedSequence]' \
--output-path aligned-seqs.qza

# import tree 
qiime tools import \
--input-path $asv_tree \
--output-path rooted-tree.qza \
--type 'Phylogeny[Rooted]'
