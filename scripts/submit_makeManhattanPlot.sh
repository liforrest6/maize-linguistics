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

module load R

results_dir='/group/jrigrp10/maize-linguistics/results'

results_files=/group/jrigrp10/maize-linguistics/results/aztecan_gemma_output/*_GWAS_results.txt

for file in $results_files
do
	Rscript makeManhattanPlot.R $file 
done