#!/bin/bash -l

#SBATCH -D /group/jrigrp10/maize-linguistics/scripts
#SBATCH -o /home/fli21/slurm-log/FEEMS-%j.txt
#SBATCH -e /home/fli21/slurm-log/FEEMS-%j.txt
#SBATCH -J FEEMS
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu 20GB
#SBATCH -n 4
#SBATCH --array=23-26
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

## requires loading conda with FEEMS
## uses FEEMS 2.0 (build from 10/2024) with custom changes to perform downsampling

module load conda
source activate feems_e

## individual submission
# python /group/jrigrp10/maize-linguistics/scripts/run_FEEMS.py -d 1 -c 0 -l 10 -n 1

## reads config file with separate parameters
config=/group/jrigrp10/maize-linguistics/scripts/FEEMS_parameters.txt

# Extract the downsample number for the current $SLURM_ARRAY_TASK_ID
dsample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Extract the cross_validate for the current $SLURM_ARRAY_TASK_ID
cross_validate=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

# Extract the lambda number for the current $SLURM_ARRAY_TASK_ID
lamb=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

# Extract the node specific variance for the current $SLURM_ARRAY_TASK_ID
node_specific_variance=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)

# Extract whether bootstrapping for the current $SLURM_ARRAY_TASK_ID
bootstrap=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)

# Extract the rep number
rep=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)

# Extract whether we are shuffling
shuffle=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $8}' $config)



# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, the same name, and the sex of the sample
echo "This is array task ${SLURM_ARRAY_TASK_ID}, \
the downsample is ${dsample}, \
the cross_validate is ${cross_validate}, \
the lambda is ${lamb}, \
the node_specific_variance is ${node_specific_variance}, \
the bootstrap parameter is ${bootstrap}, \
the rep is ${rep}, \
the shuffle parameter is ${shuffle}"

python /group/jrigrp10/maize-linguistics/scripts/run_FEEMS.py -d ${dsample} -c ${cross_validate} -l ${lamb} -n ${node_specific_variance} -b ${bootstrap} -r ${rep} -s ${shuffle} 