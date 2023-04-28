#!/bin/bash
# Author: Adam Sorbie
# Date: 11/10/2021
# Version: 0.8.6

# default parameters
# min_overlap=20 - removed until FIGARO is updated
threads=1

# fix annoying issue where script fails if path is specified without /

#while getopts a:f:r:p:o:t:m: flag
while getopts p:o:t: flag
do
  case "${flag}" in
    # a) amplicon_length=${OPTARG};;
    # f) f_primer_len=${OPTARG};;
    # r) r_primer_len=${OPTARG};;
    p) path=${OPTARG};;
    o) out=${OPTARG};;
    # m) min_overlap=${OPTARG};;
    t) threads=${OPTARG};;
    #*) echo "usage: $0 [-a] [-f] [-r] [-p] [-o] [-m] [-t]" >&2
    *) echo "usage: $0 [-p] [-o] [-t]" >&2
      exit 1 ;;
  esac
done

if ((OPTIND == 1))
then
  echo "No options specified, exiting"
  exit
fi

outdir=$(echo ${out%/})
mkdir -p $outdir

fastqc -t $threads $path/*.fastq.gz -o $outdir

multiqc $outdir -o ${outdir}/multiqc

# FIGARO
#figaro -i $path -o $out -a $amplicon_length -f $f_primer_len  -r $r_primer_len -m $min_overlap -F illumina
