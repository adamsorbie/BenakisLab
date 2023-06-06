#!/bin/bash

#SBATCH -o ./%x.%j.%N.out 
#SBATCH -D ./ 
#SBATCH -J "hammond_qc"
#SBATCH --get-user-env
#SBATCH --clusters=mpp3
#SBATCH --partition=mpp3_batch
#SBATCH --mem=80gb
#SBATCH --cpus-per-task=128
#SBATCH --mail-type=end
#SBATCH --mail-user=adam.sorbie@med.uni-muenchen.de
#SBATCH --export=NONE
#SBATCH --time=2-00:00:00

export OMP_NUM_THREADS=128
module load slurm_setup
 
# conda env
# coolMuc2
#module load anaconda3/
#eval "$(conda shell.bash hook)"
#conda activate biobakery3 

# coolMuc3
module load python/3.6_intel
eval "$(conda shell.bash hook)"
source activate biobakery3

# edit outputdir name
inputdir="/dss/dssfs02/lwp-dss-0001/pr63la/pr63la-dss-0000/ra52noz2/shotgun_data/Hammond_et_al_2022/fastq"
contam_db="/dss/dssfs02/lwp-dss-0001/pr63la/pr63la-dss-0000/ra52noz2/shotgun_data/custom_db/human_contam"
atria_exe="/dss/dsshome1/lxc05/ra52noz2/atria-3.2.2-linux/bin/atria"
# options
POLY_X=true
REMOVE_INTERMEDIATE=true

main(){
  cd $inputdir || exit
  mkdir -p merged fastqc_start fastqc_end trimmed clean 
  # cat PE reads into single file
  #cat_reads 
  # run fastqc on merged reads
  run_fastqc_start
  # run trimming and filtering
  run_trim_galore
  # remove sequences aligning to host genome
  remove_host
  # rerun fastqc after trimming and host removal
  run_fastqc_end
  # run multiqc to combined fastqc reports
  run_multiqc
  # run seqkit to get Gbp and number of seqs per read
  run_seqkit
  # clean up directories, remove unecessary files
  cleanup	
}

cat_reads(){
  # might be better to cat reads on the fly - could update this to just cat single file 
  echo "Merging reads: $(date +"%d/%m/%Y %H:%M")"
  samples=$(ls *_R1_001.fastq.gz | cut -f1 -d"_")
  for i in $samples;
  do 
    R1=${i}_S1_L001_R1_001.fastq.gz
    R2=${i}_S1_L001_R2_001.fastq.gz 
    zcat $R1 $R2 > merged/${i}_merged.fastq
  done 
  echo "Reads merged: $(date +"%d/%dm/%Y %H:%M")"
}

run_fastqc_start(){
  fastqc -t 128 merged/*.fastq -o fastqc_start
}

run_trim_galore(){
  echo "Starting trimming: $(date +"%d/%m/%Y %H:%M")"
  # loop through fastq files - need to change this or add rename script, _1 file naming pattern not most common 
  samples=$(ls merged/*.fastq | cut -f1 -d"_")
  for i in $samples;
  do 
    read=${i}_merged.fastq
    echo "trimming: ${read}"
    if [ "$POLY_X" = true ] ; then
      $atria_exe -r $read -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --length-range 75:200 -q 20 -n 2 -g NO --polyG -t 128 --force --output-dir trimmed >> trimming.log 2>&1 
      # delete untrimmed read
      rm -rf $read
    else 
      $atria_exe -r $read -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --length-range 75:200 -q 20 -n 2 -g NO -t 128 --force --output-dir trimmed >> trimming.log 2>&1 
      # delete untrimmed read
      rm -rf $read
    fi 
  done
  echo "Trimming completed: $(date +"%d/%m/%Y %H:%M")"
}

remove_host(){
  echo "Starting host removal: $(date +"%d/%m/%Y %H:%M")"
  samples=$(ls trimmed/*.fastq | cut -f1 -d".")
  for i in $samples;
  do 
    read=${i}.atria.fastq
    sample_name=$(basename $i)
    echo "decontaminating: ${read}"
    bowtie2 -x $contam_db -U $read -p 128 --very-sensitive --un clean/${sample_name}_clean.fastq 
    # delete trimmed contaminated read
    rm -rf $read 
  done
  echo "Host removal completed: $(date +"%d/%m/%Y %H:%M")"
}

run_fastqc_end(){
  fastqc -t 128 clean/*.fastq -o fastqc_end
}

run_multiqc(){
  multiqc fastqc_start -o fastqc_start
  multiqc fastqc_end -o fastqc_end
}

run_seqkit(){
  seqkit stats -j 128 -a merged/*.fastq -b -T > merged/sequence_stats_start.tsv
  seqkit stats -j 128 -a trimmed/*.fastq -b -T > trimmed/sequence_stats_trim.tsv
  seqkit stats -j 128 -a clean/*.fastq -b -T > clean/sequence_stats_end.tsv
  sequence_stats.py -r merged/sequence_stats_start.tsv \
    -t trimmed/sequence_stats_trim.tsv \
    -c clean/sequence_stats_end.tsv \
    -o sequence_stats.tsv
}
 
cleanup(){ 
  mkdir -p log reports
  # move log files to folder 
  mv *.log log
  # move sequence stats to reports folder
  mv sequence_stats.tsv reports
  if [ "$REMOVE_INTERMEDIATE" = true ] ; then
    # remove merged and trimmed files to save space
    rm -rf merged trimmed
    # move multiqc reports and delete fastqc folders
    mv fastqc_start/multiqc_report.html reports/multiqc_start.html
    mv fastqc_end/multiqc_report.html reports/multiqc_end.html
    rm -rf fastqc_start fastqc_end
  else 
    # move multiqc reports and delete fastqc folders
    cp fastqc_start/multiqc_report.html reports/multiqc_start.html
    cp fastqc_end/multiqc_report.html reports/multiqc_end.html
  fi

}

echo "Starting processing: $(date +"%d/%m/%Y %H:%M")"

main
cat reports/sequence_stats.tsv
echo "Processing complete: $(date +"%d/%m/%Y %H:%M")"
