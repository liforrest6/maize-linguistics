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

## choose levenshtein, branch, RF

Rscript all-models.R RF