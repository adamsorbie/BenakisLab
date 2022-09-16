#!/bin/bash
while getopts i:p:o: flag
 do
   case "${flag}" in
     i) input=${OPTARG};;
     p) patterns=${OPTARG};;
     o) output=${OPTARG};;
     *) echo "usage: $0 [-i] [-p] [-o]" >&2
        exit 1 ;;
   esac
done
perl -pe '/^>/ ? print "\n" : chomp' $input | grep -f $patterns -w -A1 | grep -v -- "^--$" > $output