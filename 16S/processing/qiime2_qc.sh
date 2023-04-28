#!/bin/bash 

# Author: Adam Sorbie
# Date: 14/07/22
# Version: 0.8.5 

# This script assumes you have installed the latest version of QIIME2 and activated your conda environment 
while getopts p:o:g:G:l:t: flag
 do
   case "${flag}" in
     p) path=${OPTARG};;
     o) output=${OPTARG};;
     g) f_primer=${OPTARG};;
     G) r_primer=${OPTARG};;
     l) min_len=${OPTARG};; 
     t) threads=${OPTARG};;
     *) echo "usage: $0 [-p] [-o] [-g] [-G] [-l] [-t]" >&2
        exit 1 ;;
   esac
done

# check qiime2 environment activated and exit if not 
env=$CONDA_DEFAULT_ENV
if [[ $env != *"qiime2"* ]]; then
  echo "Please activate qiime2 environment"
  exit
fi

mkdir -p $output

# create manifest file if required 
ls ${path}/*R1_001.fastq.gz | xargs -n 1 basename | cut -f1 -d"_" > q2_manifest_sn
ls ${path}/*R1_001.fastq.gz > q2_manifest_ffp
ls ${path}/*R2_001.fastq.gz > q2_manifest_rfp
paste q2_manifest_sn q2_manifest_ffp q2_manifest_rfp > q2_manifest_tmp.tsv
rm -rf q2_manifest_sn q2_manifest_ffp q2_manifest_rfp
python create_manifest.py q2_manifest_tmp.tsv ${output}/q2_manifest
rm -rf q2_manifest_tmp.tsv 

cd $output || exit 
 
# import fastq files 
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-path ${output}/q2_manifest \
  --output-path ${output}/demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2


# run cutadapt 
qiime cutadapt trim-paired \
   --i-demultiplexed-sequences demux-paired-end.qza \
   --p-front-f $f_primer \
   --p-front-r $r_primer \
   --p-minimum-length $min_len \
   --o-trimmed-sequences demux-paired-end-trimmed.qza \
   --p-cores $threads \
   --verbose


qiime demux summarize \
  --i-data ${output}/demux-paired-end-trimmed.qza \
  --o-visualization ${output}/demux-paired-end-trimmed.qzv

mkdir -p qc_out 
mv demux-paired-end-trimmed.qzv qc_out 
