#!/bin/bash

#SBATCH -o ./%x.%j.%N.out 
#SBATCH -D ./ 
#SBATCH -J "ahrtax"
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --partition=cm4_inter_large_mem
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --mail-type=end
#SBATCH --mail-user=adam.sorbie@med.uni-muenchen.de
#SBATCH --export=NONE
#SBATCH --time=01:00:00


module load slurm_setup 
export OMP_NUM_THREADS=32

module load anaconda3/
eval "$(conda shell.bash hook)"
conda activate biobakery3

# options
HIGH_SENSITIVITY=true

study="ahr_project_2023"
db="/dss/dssfs02/lwp-dss-0001/pr63la/pr63la-dss-0000/ra52noz2/shotgun_data/metaphlan_db" 
data_dir="/dss/dssfs02/lwp-dss-0001/pr63la/pr63la-dss-0000/ra52noz2/shotgun_data/AHR_project/fastq/clean"
index="mpa_vOct22_CHOCOPhlAnSGB_202212"
analysis_folder="/dss/dsshome1/lxc05/ra52noz2/shotgun_metagenomics/analysis"

main(){

    run_metaphlan
    merge_tables
	
}


run_metaphlan(){

echo "Running" $(metaphlan --v)

# make study and metaphlan folders
mkdir -p ${analysis_folder}/${study}/metaphlan/{profiles,merged_table} 

inputdatalist=$(ls -d ${data_dir}/*.fastq | awk '{print $NF}')

for read in ${inputdatalist};
do
sample=$(basename ${read} | cut -f1 -d"_")
echo "Running MetaPhlAn on sample ${sample}"
if [ "$HIGH_SENSITIVITY" = true ] ; then
    metaphlan $read --input_type fastq --min_mapq_val 2 --stat_q 0.05 --perc_nonzero 0.20 --bowtie2out $db/${sample}.bowtie2.bz2 --force --index $index --bowtie2db $db -o ${analysis_folder}/${study}/metaphlan/profiles/${sample}_profile.txt --nproc 32 &> ${analysis_folder}/${study}/${sample}.log  
  else 
    metaphlan $read --input_type fastq --bowtie2out $db/${sample}.bowtie2.bz2 --force --index $index --bowtie2db $db -o ${analysis_folder}/${study}/metaphlan/profiles/${sample}_profile.txt --nproc 32 &> ${analysis_folder}/${study}/${sample}.log 
  fi

done
 
}

merge_tables(){

echo "Merging MetaPhlAn output"
merge_metaphlan_tables.py ${analysis_folder}/${study}/metaphlan/profiles/*.txt > ${analysis_folder}/${study}/metaphlan/merged_table/merged_table.txt

}

echo "Starting analysis: $(date +"%m/%d/%Y %H:%M")"
main
echo "Analysis complete: $(date +"%m/%d/%Y %H:%M")"
