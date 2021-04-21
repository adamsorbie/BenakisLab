#!/bin/bash
# Author: Adam Sorbie 
# Date: 15/04/21
# Version 0.6.0

# $1 = delim
# $2 = last field
# $3 = fwd pattern e.g. R1.fastq
# $4 = rev pattern e.g. R1.fastq
# rename forward reads
for i in $2;
do
    name=$(echo $i | cut -f1-$2 -d "$1")
	mv $i ${name}"-S1_L001_R1_001.fastq"
done
# rename reverse reads
for i in *R.fastq;
do
    name=$(echo $i | cut -f1-$2 -d "$1")
	mv $i ${name}"-S1_L001_R2_001.fastq"
done
