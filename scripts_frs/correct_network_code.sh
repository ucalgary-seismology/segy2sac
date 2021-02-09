#!/bin/bash

#SBATCH --partition=geo,glados12,glados16
#SBATCH --time=1-00:00:00
#SBATCH --mem=3G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --array=1-305%40
#SBATCH --output="outslurm/slurm-%A_%a.out"

echo "Starting task $SLURM_ARRAY_TASK_ID"
path=$(sed -n "${SLURM_ARRAY_TASK_ID}p" paths.list)
#path="/home/gilbert_lab/cami_frs/borehole_data/sac_daily_500Hz/20200108"

echo $path
curdir=`pwd`
cd $path

for f in `ls *.sac`
do 
	network=`echo $f | awk -F. '{print $1}'`
	station=`echo $f | awk -F. '{print $2}'`
	echo "Current network & station: $network $station"
#	newstation=`echo $f | awk -F. '{printf "%s%04d", substr($2,1,1), substr($2,2)}'`
        newstation=`echo $f | awk -F. '{printf "BH%03d", substr($2,2)}'`
	echo "New network station: 8O $newstation"
	
	if [[ "$network" == "BH" ]]
	then
		# First change network code and station code
		printf "r $f\nch knetwk 8O kstnm $newstation\nwh\nq\n" | sac

		# Then rename file
		newf=`echo $f | sed "s/${network}/8O/" | sed "s/${station}/${newstation}/"`
		echo $f $newf
		mv $f $newf
	fi
done

cd $curdir
