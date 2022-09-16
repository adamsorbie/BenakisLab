#!/bin/bash 

# normalize qiime artifiact with minimum sum scaling   
python normalize_mss.py "$@" 

# re-import
biom convert -i feature-table-norm.tsv -o feature-table-norm.biom --table-type="OTU table" --to-hdf5

qiime tools import \
  --input-path feature-table-norm.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature-table-norm.qza

  # cleanup 
  rm -rf feature-table-norm.tsv feature-table-norm.biom 