#!/bin/bash

# Author: Adam Sorbie
# Date: 20/12/21
# Version: 0.1.1 

while getopts p: flag
do
  case "${flag}" in
    p) path=${OPTARG};; 
    *) echo "usage: $0 [-p]" >&2
      exit 1 ;;
  esac
done

cd $path 

filenames=$(ls *R1_001.fastq.gz | cut -f1 -d"_")

mkdir phred33 

for i in $filenames;
do 
  R1=${i}_S1_L001_R1_001.fastq.gz
  R2=${i}_S1_L001_R2_001.fastq.gz
  out1=phred33/${R1}
  out2=phred33/${R2} 
  reformat.sh in=$R1 in2=$R2 out=$out1 out2=$out2 qin=64 qout=33
done




