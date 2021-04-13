#!/bin/bash 

## required options - 
# -a
# -A
# -m
# -M
# -o 

while getopts a:A:m:M:o: flag
do
  case "${flag}" in
    a) f_primer=${OPTARG};;
    A) r_primer=${OPTARG};;
    m) min_lenf=${OPTARG};;
    M) min_lenr=${OPTARG};;
    o) output=${OPTARG};;
  esac
done

if ((OPTIND == 1))
then
  echo "No options specified, exiting"
  exit 
fi 



conda activate bioinfo

mkdir trimmed_primer

 for sample in $(ls *.fastq.gz);
 do;
    echo "Trimming sample: $sample"
    cutadapt -a $f_primer \
    -A $r_primer \
    -m $m -M $M --discard-untrimmed \
    -o ${sample}_R1_trimmed.fastq.gz -p ${sample}_R2_trimmed.fastq.gz \
     ${sample}_sub_R1.fq ${sample}_sub_R2.fq \
     >> cutadapt_primer_trimming_stats.txt 2>&1
done
