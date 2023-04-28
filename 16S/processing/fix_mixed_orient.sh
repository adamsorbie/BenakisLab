#!/bin/bash 
# Author: Adam Sorbie 
# Date:22/07/21
# Version: 0.1.0

while getopts g:G:m:M:p: flag
do
  case "${flag}" in
    g) f_primer=${OPTARG};;
    G) r_primer=${OPTARG};;
    m) min_len=${OPTARG};;
    M) max_len=${OPTARG};;
    p) path=${OPTARG};;
    *) echo "usage: $0 [-g] [-G] [-m] [-M] [-p]" >&2
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

mkdir -p cat_dupes 

ls *_R1_001.fastq.gz | cut -f1 -d"_" > sample_names 

while read i; do
  R1=${i}_S1_L001_R1_001.fastq.gz
  R2=${i}_S1_L001_R2_001.fastq.gz
  cat $R1 $R2 > cat_dupes/${i}_R1R2.fastq.gz 
  cat $R2 $R1 > cat_dupes/${i}_R2R1.fastq.gz
done < sample_names

eval "$(conda shell.bash hook)"
conda activate bioinfo

cd cat_dupes 

mkdir -p trimmed_primer

# this currently loops through all files rather than pairs of files 
for sample in $(ls *R1R2.fastq.gz | cut -f1 -d_);
do
    echo "Trimming sample: $sample"
    cutadapt -g $f_primer -G $r_primer -m $min_len -M $max_len --discard-untrimmed -j 0 -o ${sample}_trimmed_primer_R1_001.fastq.gz -p ${sample}_trimmed_primer_R2_001.fastq.gz ${sample}_R1R2.fastq.gz ${sample}_R2R1.fastq.gz >> trimmed_primer/cutadapt_primer_trimming_stats.txt 2>&1
done

mv *trimmed_primer*.fastq.gz trimmed_primer
