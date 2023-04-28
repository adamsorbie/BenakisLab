#!/bin/bash
# $1 input fasta
# $2 patterns file
# $3 output fasta
dos2unix $2
perl -pe '/^>/ ? print "\n" : chomp' $1 | grep -f $2 -w -A1 | grep -v -- "^--$" > $3
