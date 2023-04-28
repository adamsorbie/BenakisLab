#!/bin/bash 

# $1 = fastq directory
cd $1 
echo $PWD

for i in *.fastq.gz; 
do
  zcat $i | awk '{if(NR%4==2) print length($1)}' >  ${i}.readslength.txt
done

mkdir read_length_distribution
mv *.readslength.txt read_length_distribution
