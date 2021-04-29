#!/bin/bash
# Author: Adam Sorbie
# Date: 26/04/2021
# Version: 0.6.0
 
# default 
min_overlap=20

while getopts a:f:r:p:o:m: flag
do
  case "${flag}" in
    a) amplicon_length=${OPTARG};;
    f) f_primer_len=${OPTARG};;
    r) r_primer_len=${OPTARG};;
    p) path=${OPTARG};;
    o) out=${OPTARG};;
    m) min_overlap=${OPTARG};;
    *) echo "usage: $0 [-a] [-f] [-r] [-p] [-o]" >&2
       exit 1 ;;
  esac
done

if ((OPTIND == 1))
then
  echo "No options specified, exiting"
  exit
fi


# activate conda env with fastqc installed
eval "$(conda shell.bash hook)"
conda activate bioinfo

fastqc -t 8 $path/*.fastq.gz -o $out

multiqc $out -o ${out}/multiqc
# FIGARO
figaro -i $path -o $out -a $amplicon_length -f $f_primer_len  -r $r_primer_len -m $min_overlap -F illumina 

