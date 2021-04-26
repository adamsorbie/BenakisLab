#!/bin/bash
# Author: Adam Sorbie
# Date: 26/04/2021
# Version: 0.5.0
# $1 raw data directory 
# $2 out directory
# need getopts here path, output, amplicon length, f and r primer length

# Input director contains fastq
raw_data=$1

if [[ $# -eq 0 ]] ; then
    echo 'No arguments provided'
    exit 0
fi

fastqc -t 8 $raw_data/*.fastq.gz 
cd $raw_data
mkdir -p $2 
mv *fastqc* $2

multiqc $2 

# FIGARO
figaro -i  -o /path/to/output/files -a [amplicon length] \
    -f [forward primer length] -r [reverse primer length]
