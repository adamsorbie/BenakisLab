#!/bin/bash 

# Author: Adam Sorbie 
# Date: 15/04/20
# Version: 0.5.0

while getopts a:A:m:M: flag
do
  case "${flag}" in
    a) f_primer=${OPTARG};;
    A) r_primer=${OPTARG};;
    m) min_lenf=${OPTARG};;
    M) min_lenr=${OPTARG};;
    *) echo "usage: $0 [-a] [-A] [-m] [-M]" >&2
       exit 1 ;;
  esac
done

if ((OPTIND == 1))
then
  echo "No options specified, exiting"
  exit 
fi 

conda activate bioinfo

mkdir trimmed_primer

for sample in $(ls *.fastq.gz | cut -f1-3 -d"_");
do
    echo "Trimming sample: $sample"
    cutadapt -a $f_primer \
    -A $r_primer \
    -m $min_lenf -M $min_lenr --discard-untrimmed \
    -o ${sample}_R1_trimmed_primer.fastq.gz -p ${sample}_R2_trimmed_primer.fastq.gz \
     ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz \
     >> cutadapt_primer_trimming_stats.txt 2>&1
done

mkdir trimmed_primer 
mv *trimmed_primer.fastq.gz trimmed_primer
