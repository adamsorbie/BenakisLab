#!/bin/bash 

# Author: Adam Sorbie
# Date: 14/09/21
# Version: 0.0.1 

# This script assumes you have installed the latest version of QIIME2 and activated your conda environment 

# check qiime2 environment activated and exit if not 
env=$CONDA_DEFAULT_ENV
if [[ $env != *"qiime2"* ]]; then
  echo "Please activate qiime2 environment"
  exit
fi

mkdir qiime2_pipeline && cd "$_"

# import fastq files 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $path \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

qiime demux summarize \
    --i-data demux-paired-end.qza \
    --o-visualization ./demux-paired-end.qzv


