#!/bin/bash -l

#SBATCH -D /group/jrigrp10/maize-linguistics/scripts
#SBATCH -o /home/fli21/slurm-log/linguistics-%j.txt
#SBATCH -e /home/fli21/slurm-log/linguistics-%j.txt
#SBATCH -J linguistics
#SBATCH -t 4:00:00
#SBATCH --mem 10GB
#SBATCH -n 8
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

module load gemma

results_dir='/group/jrigrp10/maize-linguistics/results/'
data_dir='/group/jrigrp10/maize-linguistics/data/'

language='otomanguean'

if [ $language == 'mayan' ]; then
	groupings=("Chicomuceltec"
	"Chol"
	"Chuj"
	"Yucatec"
	"Tzutujil"
	"Tzotzil");
	num_groupings=5
elif [ $language == 'aztecan' ]; then
	groupings=("Tepecano"
	"Pochutec"
	"Nahuatl"
	"Cora"
	"Pipil"
	"Tahue")
	num_groupings=5
elif [ $language == 'otomanguean' ]; then
	groupings=("Mixtec"	
	"Matlatzinca"	
	"Zapotec"	
	"Mangue"	
	"Chinantecan"	
	"Popolocan"	
	"Otomi"	
	"Pame")
	num_groupings=6
fi


for i in $(seq 0 $num_groupings);
do 

	# gemma -g ${data_dir}/${language}/${language}.dosage.gemma.vcf \
	# -p ${data_dir}/${language}/${language}_phenotype.txt  \
	# -k ${data_dir}/${language}/Kinship_filtered_cleaned.csv \
	# -lmm 4 \
	# -n $((i+1))\
	# -outdir /group/jrigrp10/maize-linguistics/results/${language}_gemma_output \
	# -o ${groupings[$i]}_results

	cut -f2- /group/jrigrp10/maize-linguistics/results/${language}_gemma_output/${groupings[$i]}_results.assoc.txt | \
	awk 'NR==FNR{a[$2]=$1; b[$2]=$3; next}{OFS="\t"}{print b[$1], a[$1], $0}' /group/jrigrp10/maize-linguistics/data/${language}/${language}_chromosome.txt - | \
	# awk -F'\t' '{print $1"_"$2 FS $0}' -  | \
	cut -f1-3,16 - | \
	sed '1c\SNP\tCHR\tBP\tP' > /group/jrigrp10/maize-linguistics/results/${language}_gemma_output/${groupings[$i]}_GWAS_results.txt

done

	# cut -f2- /group/jrigrp10/maize-linguistics/results/mayan_gemma_output/Yucatec_results.assoc.txt | \
	# awk 'NR==FNR{a[$2]=$1; b[$2]=$3; next}{OFS="\t"}{print b[$1], a[$1], $0}' /group/jrigrp10/maize-linguistics/data/mayan/mayan_chromosome.txt - | \
	# # awk -F'\t' '{print $1"_"$2 FS $0}' -  | \
	# cut -f1-3,16 - | \
	# sed '1c\SNP\tCHR\tBP\tP' | head
