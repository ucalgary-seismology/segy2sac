#!/bin/bash

#SBATCH -p geo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=1-27
#SBATCH --output="outslurm/slurm-%A_%a.out"

output_dir="/home/gilbert_lab/cami_frs/borehole_data/sac_hourly_nez_500Hz"
input_dir="/home/gilbert_lab/cami_frs/borehole_data/sac_daily_nez_500Hz"

day=$(sed -n "${SLURM_ARRAY_TASK_ID}p" folders.list)
echo $day

output_path="%n.%s.%l.%c.%(wmin)s.%(wmax)s.sac"
#jackseis ${input_dir}/${day} --tinc=3600 --downsample=125 --output-dir=${output_dir}/${day} --format=sac --output=${output_path} --output-format=sac
jackseis ${input_dir}/${day} --force --tinc=3600 --downsample=125 --output-dir=${output_dir}/${day} --format=sac --output=${output_path} --output-format=sac
