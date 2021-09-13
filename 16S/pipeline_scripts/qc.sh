#!/bin/bash
# Author: Adam Sorbie
# Date: 30/04/2021
# Version: 0.7.5
 
# default 
min_overlap=20

# fix annoying issue where script fails if path is specified without /

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
outdir=$(echo ${out%/}) 

fastqc -t 8 $path/*.fastq.gz -o $outdir

multiqc $outdir -o ${outdir}/multiqc
# FIGARO
figaro -i $path -o $out -a $amplicon_length -f $f_primer_len  -r $r_primer_len -m $min_overlap -F illumina 

