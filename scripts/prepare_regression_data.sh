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

## author: Forrest Li
## script to prepare file formats for linguistic regression

module load plink
# module load vcftools
module load bcftools

dataset='Haynie'

results_dir=/group/jrigrp10/maize-linguistics/results/${dataset}
data_dir='/group/jrigrp10/maize-linguistics/data/'

language='mayan'

## prepare sample list for plink format ingestion
cat ${results_dir}/${language}_${dataset}_accessions.txt |
uniq | 
sort |
awk -v FS='\t' -v OFS='\t' 'BEGIN {OFS = FS} {print $0 FS $1}' |
uniq > ${data_dir}/${language}/${language}_${dataset}_accessions_plink_format.txt

## plink to do LD pruning to r2 = 0.5
plink --vcf /group/jrigrp10/smambakk/language/zeaGBS20161020ZeaAlbertoTaxaAGPv4_forBeagle.imputed.vcf.vcf \
--keep ${data_dir}/${language}/${language}_${dataset}_accessions_plink_format.txt \
--allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--maf 0.01 \
--export vcf \
--out ${data_dir}/${language}/${language}.${dataset}.regression


# ## get dosage file format for regressions
bcftools +dosage ${data_dir}/${language}/${language}.${dataset}.regression.vcf -i ID=@${data_dir}/${language}/${language}.${dataset}.regression.prune.in > ${data_dir}/${language}/${language}.${dataset}.LD-pruned.dosage.vcf
cut -f5- ${data_dir}/${language}/${language}.${dataset}.LD-pruned.dosage.vcf | sed -E 's/([0-9]+)\.0\b/\1/g' > ${data_dir}/${language}/${language}.${dataset}.LD-pruned.dosage.regression.vcf #Turn your VCF into dosage


# bcftools +dosage /group/jrigrp10/maize-linguistics/data//mayan/mayan.NL.regression.vcf -i ID=@/group/jrigrp10/maize-linguistics/data//mayan/mayan.NL.regression.prune.in
# ## remove intermediary files
rm ${data_dir}/${language}/${language}.${dataset}.LD-pruned.dosage.vcf
rm ${data_dir}/${language}/${language}.${dataset}.regression.vcf
rm ${data_dir}/${language}/${language}.${dataset}.regression.nosex
rm ${data_dir}/${language}/${language}.${dataset}.regression.log

# grep “^#CHROM” MAF0.01_21km_MaizeGBS.vcf > LDprune_MAF0.01_21km_MaizeGBS.tsv
# grep -w -f maize.prune.in MAF0.01_21km_MaizeGBS.vcf  >> LDprune_MAF0.01_21km_MaizeGBS.tsv
# #then change genotype scores to be 0, 1, 2 (numeric)
# sed -i ‘s/0\/0/0/g’ LDprune_MAF0.01_21km_MaizeGBS.tsv
# sed -i ‘s/0\/1/1/g’ LDprune_MAF0.01_21km_MaizeGBS.tsv
# sed -i ‘s/1\/0/1/g’ LDprune_MAF0.01_21km_MaizeGBS.tsv
# sed -i ‘s/1\/1/2/g’ LDprune_MAF0.01_21km_MaizeGBS.tsv
# sed -i ‘s/\.\/\./NA/g’ LDprune_MAF0.01_21km_MaizeGBS.tsv