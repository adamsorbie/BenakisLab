#!/bin/bash

#SBATCH -o ./%x.%j.%N.out 
#SBATCH -D ./ 
#SBATCH -J "ahr_func"
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --partition=cm4_inter_large_mem
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128gb
#SBATCH --mail-type=end
#SBATCH --mail-user=adam.sorbie@med.uni-muenchen.de
#SBATCH --export=NONE
#SBATCH --time=72:00:00

export OMP_NUM_THREADS=48
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate biobakery3

study="ahr_project_2023"
data_dir="/dss/dssfs02/lwp-dss-0001/pr63la/pr63la-dss-0000/ra52noz2/shotgun_data/AHR_project/fastq/clean"
analysis_folder="/dss/dsshome1/lxc05/ra52noz2/shotgun_metagenomics/analysis"
metaphlan_profile="/dss/dsshome1/lxc05/ra52noz2/shotgun_metagenomics/analysis/ahr_project_2023/metaphlan/merged_table/merged_table.txt"

main(){
    
	run_humann
	normalize
	merge_tables
  regroup_tables
  generate_summed_tables 
}


run_humann(){


echo "Running" $(humann --version)

# make metaphlan folder
mkdir -p ${analysis_folder}/${study}/humann/{profiles,merged_table} 

cd ${analysis_folder}/${study}/humann || exit
filelist=$(ls -d ${data_dir}/*.fastq | awk '{print $NF}')

for read in ${filelist};
do
sample=$(basename ${read} | cut -f1 -d"_") 
# run humann 
humann --input $read --output ${analysis_folder}/${study}/humann/profiles/ --taxonomic-profile $metaphlan_profile --threads 48 --remove-temp-output --o-log ${analysis_folder}/${study}/humann/${fname}.log
done

}

normalize(){

genefamilies=$(ls -d ${analysis_folder}/${study}/humann/profiles/*genefamilies* | awk '{print $NF}')

for i in $genefamilies;
do 
  sample=$(basename ${i} | cut -f1 -d"_") 
  humann_renorm_table --input $i --output ${analysis_folder}/${study}/humann/profiles/${sample}_genefamilies_relab.tsv --units relab
done

pathabundance=$(ls -d ${analysis_folder}/${study}/humann/profiles/*pathabundance* | awk '{print $NF}')

for i in $pathabundance;
do 
  sample=$(basename ${i} | cut -f1 -d"_") 
  humann_renorm_table --input $i --output ${analysis_folder}/${study}/humann/profiles/${sample}_pathabundance_relab.tsv --units relab
done



}

merge_tables(){

echo "Merging HUMAnN output"
humann_join_tables --input ${analysis_folder}/${study}/humann/profiles/ --output ${analysis_folder}/${study}/humann/merged_table/humann_genefamilies_relab.tsv --file_name genefamilies_relab
humann_join_tables --input ${analysis_folder}/${study}/humann/profiles/ --output ${analysis_folder}/${study}/humann/merged_table/humann_pathcoverage_relab.tsv --file_name pathcoverage
humann_join_tables --input ${analysis_folder}/${study}/humann/profiles/ --output ${analysis_folder}/${study}//humann/merged_table/humann_pathabundance_relab.tsv --file_name pathabundance_relab

}

regroup_tables(){
  echo "Regrouping HUMAnN output"
  humann_regroup_table --input ${analysis_folder}/${study}/humann/merged_table/humann_genefamilies_relab.tsv --groups uniref90_ko  --output ${analysis_folder}/${study}/humann/merged_table/humann_genefamilies_ko.tsv
  humann_regroup_table --input ${analysis_folder}/${study}/humann/merged_table/humann_genefamilies_relab.tsv --groups uniref90_level4ec --output ${analysis_folder}/${study}/humann/merged_table/humann_genefamilies_ec.tsv
}

generate_summed_tables(){
 echo "Creating unstratified tables"
 humann_split_stratified_table -i ${analysis_folder}/${study}/humann/merged_table/humann_genefamilies_ko.tsv -o ${analysis_folder}/${study}/humann/merged_table/unstrat
 humann_split_stratified_table -i ${analysis_folder}/${study}/humann/merged_table/humann_genefamilies_ec.tsv -o ${analysis_folder}/${study}/humann/merged_table/unstrat
 # remove stratified 
 rm -rf ${analysis_folder}/${study}/humann/merged_table/unstrat/*_stratified.tsv 
}

main
echo "Analysis complete"
