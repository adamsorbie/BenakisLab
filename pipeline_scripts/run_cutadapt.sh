#!/bin/bash 

# Author: Adam Sorbie 
# Date: 15/04/20
# Version: 0.8.2

while getopts a:A:m:M:p: flag
do
  case "${flag}" in
    a) f_primer=${OPTARG};;
    A) r_primer=${OPTARG};;
    m) min_len=${OPTARG};;
    M) max_len=${OPTARG};;
    p) path=${OPTARG};;
    *) echo "usage: $0 [-a] [-A] [-m] [-M] [-p]" >&2
       exit 1 ;;
  esac
done

if ((OPTIND == 1))
then
  echo "No options specified, exiting"
  exit 
fi 

cd $path || 

conda activate bioinfo

mkdir -p trimmed_primer

for sample in $(ls *.fastq.gz | cut -f1-3 -d"_");
do
    echo "Trimming sample: $sample"
    cutadapt -a $f_primer \
    -A $r_primer \
    -m $min_len -M $max_len \
    -o ${sample}_trimmed_primer_R1_001.fastq.gz -p ${sample}_trimmed_primer_R2_001.fastq.gz \
     ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz \
     >> trimmed_primer/cutadapt_primer_trimming_stats.txt 2>&1
done

mv *trimmed_primer*.fastq.gz trimmed_primer
