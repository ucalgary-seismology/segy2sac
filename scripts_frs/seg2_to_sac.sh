#!/bin/bash

#SBATCH --array=1-330
#SBATCH --nodes=1
#SBATCH --partition=geo,glados12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=1-00:00:00
#SBATCH --output="outslurm/slurm-%A_%a.out"

echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
folder=`sed -n "${SLURM_ARRAY_TASK_ID}p" alldays.list`
echo $folder

python process_seg2_to_sac.py $folder
#python write_sac_to_NEZ.py $folder
