#!/bin/bash

#SBATCH --partition=geo,glados12
#SBATCH --time=10:00:00
#SBATCH --array=1-49
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output="outslurm/slurm-%A_%a.out"

echo "Starting task $SLURM_ARRAY_TASK_ID"
day=$(sed -n "${SLURM_ARRAY_TASK_ID}p" dates_rotate.list)

python write_sac_to_NEZ.py $day
