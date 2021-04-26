#!/bin/bash 

for i in $@;
do
   sed -i 's/Zotu/ASV_/g' $i
done
