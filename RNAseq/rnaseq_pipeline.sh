#!/bin/bash 

#SBATCH -J str_sh_rnaseq
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std 
#SBATCH --mail-type=end
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=adam.sorbie@med.uni-muenchen.de
#SBATCH --export=NONE
#SBATCH --time=8:00:00

export OMP_NUM_THREADS=16 

source /etc/profile.d/modules.sh

module load python/3.6_intel 
source activate rnaseq 

datadir="/dss/dssfs02/lwp-dss-0001/pr63la/pr63la-dss-0000/ra52noz2/rawdata"
mkdir -p qc/{fastqc,multiqc}

# Run fastqc 
fq_ext=".fq.gz"
fastqc -o qc/fastqc -t 16 "${datadir}"/*"${fq_ext}" 

# Build salmon index 
## check if it exists and skip if yes
salmon_index="mouse_transcriptome/mouse_index.idx"
if [ -f $salmon_index ]; then 
  echo "Index already exists, skipping index building" 
else
  ## make directory for transcriptome
  mkdir -p transcriptome && cd $_

  ## download mouse transcriptome inc ncRNA and combine into one fasta file
  wget ftp://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/cdna/*.fa.gz .
  wget ftp://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/ncrna/*fa.gz .
  zcat *.fa.gz > mouse_transcriptome.fa
  ## build salmon index file
  salmon index -t mouse_transcriptome.fa -i mouse_index.idx -p 3 -k 31
  cd ..
fi 

# Run salmon pseudoalignment
file_ext="R1_001.fq.gz"

for fn in "${datadir}"/*"${file_ext}";
do
  samp=$(basename $fn | cut -f1,2 -d"_")
  echo "Processing sample ${samp}"
  salmon quant -i transcriptome/mouse_index.idx -l A -1 ${datadir}/${samp}_S1_L001_R1_001.fq.gz -2 ${datadir}/${samp}_S1_L001_R2_001.fq.gz -p 16 -o quants/${samp}_quant
done

# run multiqc 
multiqc . -o qc/multiqc -n rnaseq_report	     
