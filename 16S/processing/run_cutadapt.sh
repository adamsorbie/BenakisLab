#!/bin/bash

# Author: Adam Sorbie
# Date: 01/06/21
# Version: 0.8.3

# default parameters
paired_end=true

while getopts g:G:m:M:p:P: flag
do
  case "${flag}" in
    g) f_primer=${OPTARG};;
    G) r_primer=${OPTARG};;
    m) min_len=${OPTARG};;
    M) max_len=${OPTARG};;
    p) path=${OPTARG};;
    P) paired_end=${OPTARG};;
    *) echo "usage: $0 [-g] [-G] [-m] [-M] [-p] [-P]" >&2
       exit 1 ;;
  esac
done


if ((OPTIND == 1))
then
  echo "No options specified, exiting"
  exit
fi

echo $path
cd $path ||

echo $PWD

eval "$(conda shell.bash hook)"
conda activate bioinfo

mkdir -p trimmed_primer

if [ "$paired_end" = true ]; then
  for sample in $(ls *.fastq.gz | cut -f1 -d"_");
  do
    echo "Trimming sample: $sample"
    cutadapt -g $f_primer \
    -G $r_primer \
    -m $min_len -M $max_len \
    -j 0 \
    -o ${sample}_trimmed_primer_R1_001.fastq.gz -p ${sample}_trimmed_primer_R2_001.fastq.gz \
     ${sample}_S*_L001_R1_001.fastq.gz ${sample}_S*_L001_R2_001.fastq.gz \
     >> trimmed_primer/cutadapt_primer_trimming_stats.txt 2>&1
  done
elif [ "$paired_end" = false ]; then
  for sample in $(ls *.fastq.gz | cut -f1 -d"_");
  do
    echo "Trimming sample: $sample"
    cutadapt -g $f_primer \
    -m $min_len -M $max_len \
    -j 0 \
    -o ${sample}_trimmed_primer_R1_001.fastq.gz  \
     ${sample}_S*_L001_R1_001.fastq.gz  \
     >> trimmed_primer/cutadapt_primer_trimming_stats.txt 2>&1
  done
fi
mv *trimmed_primer*.fastq.gz trimmed_primer
