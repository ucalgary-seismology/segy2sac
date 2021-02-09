#!/bin/bash

#SBATCH --array=1-21
#SBATCH --nodes=1
#SBATCH --partition=geo,glados12,glados16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=12:00:00
#SBATCH --output="outslurm/slurm-%A_%a.out"

echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
folder=`sed -n "${SLURM_ARRAY_TASK_ID}p" folders_segy.list`
echo $folder

python rearrange_seg2_directories.py $folder

