#!/bin/bash

#SBATCH --array=1-5
#SBATCH --nodes=1
#SBATCH --partition=geo,glados12,glados16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --output="outslurm/slurm-%A_%a.out"

echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
folder=`sed -n "${SLURM_ARRAY_TASK_ID}p" folders_rearrange.list`
echo $folder

data_dir="/home/gilbert_lab/cami_frs/borehole_data/sac_daily_nez_500Hz/"
python correct_units.py $folder $data_dir
