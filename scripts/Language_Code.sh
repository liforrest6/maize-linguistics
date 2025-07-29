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
module load vcftools
module load bcftools
# module load R

#Rscript /group/jrigrp10/maize-linguistics/scripts/One_hot_encoding.R
#echo '0_1 done'

dataset='NL'

results_dir='/group/jrigrp10/maize-linguistics/results/${dataset}'
data_dir='/group/jrigrp10/maize-linguistics/data/'

language='mayan'

## Making phenotype file
# Extracts the first column of this file, which tells you which samples go into the phenotype file
cut -f1 ${data_dir}/${language}/0_1_${language}.txt > ${data_dir}/${language}/${language}_sample_list.txt  
# Removing first two lines for just phenotype
cut -f 2- ${data_dir}/${language}/0_1_GWAS_${language}.txt | sed '1d' > ${data_dir}/${language}/${language}_phenotype.txt 

## Making dosage file
# subset vcf file
vcftools --vcf /group/jrigrp10/smambakk/language/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--keep ${data_dir}/${language}/${language}_sample_list.txt \
--recode \
--recode-INFO-all \
--out ${data_dir}/${language}/${language}.v4.imputed

sed 1,10d ${data_dir}/${language}.v4.imputed.recode.vcf | cut -f 1-3 > ${data_dir}/${language}_chromosome.txt

bcftools +dosage ${data_dir}/${language}/${language}.v4.imputed.recode.vcf > ${data_dir}/${language}/${language}.dosage.vcf
cut -f2- ${data_dir}/${language}/${language}.dosage.vcf | sed '1d' > ${data_dir}/${language}/${language}.dosage.gemma.vcf #Turn your VCF into dosage

#Making kinship matrix
grep -wFf ${data_dir}/${language}/${language}_sample_list.txt  ${data_dir}/K_allChr.csv  > ${data_dir}/${language}/Kinship_filtered.csv #Filter for the unique IDs on columns
awk -f /group/jrigrp10/maize-linguistics/scripts/transpose.awk ${data_dir}/${language}/Kinship_filtered.csv \
| grep -wFf ${data_dir}/${language}/${language}_sample_list.txt - \
| awk -f /group/jrigrp10/maize-linguistics/scripts/transpose.awk - \
| sed '1d'| cut -d',' -f2- \
> ${data_dir}/${language}/Kinship_filtered_cleaned.csv

