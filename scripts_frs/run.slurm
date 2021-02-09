#!/bin/bash

#SBATCH --partition=geo
#SBATCH --time=7-00:00:00
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output="outslurm/slurm-%A_%a.out"
#SBATCH --array=1-27

# CSV station file path
stafile="/home/gilbert_lab/cami_frs/station_deployments_info/X16_Nov-Dec2020_data.csv"

# Output parent directory
outputdir="/home/gilbert_lab/cami_frs/hawk_data/sac_data_mps/X16_Nov-Dec2020_sac_daily_500Hz/"

# Get segy file path from file "files.list"
sgyfile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" files.list)
echo "SEGY file: ${sgyfile}"

python process_segy.py $sgyfile $stafile $outputdir


