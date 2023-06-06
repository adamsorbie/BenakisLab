#!/bin/bash

#SBATCH -o ./%x.%j.%N.out 
#SBATCH -D ./ 
#SBATCH -J "shotgun_data"
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --cpus-per-task=12
#SBATCH --mem=32gb
#SBATCH --mail-type=end
#SBATCH --mail-user=adam.sorbie@med.uni-muenchen.de
#SBATCH --export=NONE
#SBATCH --time=24:00:00

cd /dss/dssfs02/lwp-dss-0001/pr63la/pr63la-dss-0000/ra52noz2/

module load miniconda3
source activate grabseqs

cd shotgun_data
while read a; do
  prid=$a
  echo $prid
  grabseqs sra $a -r 3 -t 12 -m "metadata_${prid}.csv" -o ${prid}_data/ 
done < prjids_stroke.txt
