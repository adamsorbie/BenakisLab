#!/bin/bash 

# Author: Adam Sorbie
# Date: 27/07/22
# Version: 0.9.2 

threads=1
maxeeF=2
maxeeR=2
while getopts p:g:G:m:l:f:r:n:N:t: flag
do
  case "${flag}" in
    p) path=${OPTARG};; 
    m) metadata=${OPTARG};;
    f) trunc_f=${OPTARG};;
    r) trunc_r=${OPTARG};;
    n) maxeeF=${OPTARG};;
    N) maxeeR=${OPTARG};;
    t) threads=${OPTARG};; 
    *) echo "usage: $0 [-p] [-g] [-G] [-m] [-l] [-f] [-r] [-n] [-N] [-t]" >&2
      exit 1 ;;
  esac
done

# This script assumes you have installed the latest version of QIIME2 and activated your conda environment 

# check qiime2 environment activated and exit if not 
env=$CONDA_DEFAULT_ENV
if [[ $env != *"qiime2"* ]]; then
  echo "Please activate qiime2 environment"
  exit
fi
cp filter_low_abundant.py filter_fasta.sh ${path}
cd $path || exit


qiime dada2 denoise-paired \
   --i-demultiplexed-seqs demux-paired-end-trimmed.qza \
   --p-trunc-len-f $trunc_f \
   --p-trunc-len-r $trunc_r \
   --p-max-ee-f $maxeeF \
   --p-max-ee-r $maxeeR \
   --p-trunc-q 3 \
   --p-chimera-method 'consensus' \
   --p-n-reads-learn 100000 \
   --p-n-threads $threads \
   --o-table table.qza \
   --o-representative-sequences rep-seqs.qza \
   --o-denoising-stats denoising-stats.qza \
   --verbose

#Visualize DADA2 stats
qiime metadata tabulate \
--m-input-file denoising-stats.qza \
--o-visualization denoising-stats.qzv

## Write out abundance table and filter by 0.25% 

qiime tools export --input-path table.qza --output-path .
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv

qiime tools export --input-path rep-seqs.qza --output-path .

# Filter 
# remove constructed from biom file header
tail -n +2 feature-table.tsv > feature-table-2.tsv
python filter_low_abundant.py   
bash filter_fasta.sh -i dna-sequences.fasta -p keep-asvs.txt -o rep-seqs-filt.fasta
# remove copied scripts 
rm -rf filter_low_abundant.py filter_fasta.sh 
# Re-convert to biom 
biom convert -i feature-table-filt.tsv -o feature-table-filt.biom --table-type="OTU table" --to-hdf5

## Re-import filtered table and fasta file
qiime tools import \
  --input-path feature-table-filt.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature-table-filt.qza

qiime tools import \
  --input-path rep-seqs-filt.fasta \
  --output-path rep-seqs-filt.qza \
  --type 'FeatureData[Sequence]'


## Taxonomy assignment 
if [ ! -f silva-138-99-nb-classifier.qza ];
then 
  wget \
    -O "silva-138-99-nb-classifier.qza" \
    "https://data.qiime2.org/2021.11/common/silva-138-99-nb-classifier.qza"
fi

qiime feature-classifier classify-sklearn \
  --i-reads rep-seqs-filt.qza \
  --i-classifier silva-138-99-nb-classifier.qza \
  --p-confidence 0.5 \
  --p-reads-per-batch 5000 \
  --verbose \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-filt.qza \
  --o-visualization rep-seqs-filt.qzv

## Tree generation 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-filt.qza \
  --o-alignment rep-seqs-filt-aln.qza \
  --o-masked-alignment rep-seqs-filt-aln-mask.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

## feature table summary 
qiime feature-table summarize \
--i-table feature-table-filt.qza \
--m-sample-metadata-file $metadata \
--o-visualization feature-table-filt.qzv

## TO-DO add tidy up - place files into folders
mkdir -p q2_out intermediate logs

mv table.qza demux-paired-end.qza demux-paired-end-trimmed.qza keep-asvs.txt *.tsv *.fasta *.biom
mv denoising* logs 