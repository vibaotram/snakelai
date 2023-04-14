#!/bin/bash
#SBATCH --job-name average_elai
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=long
#SBATCH --error /shared/projects/elai_most/data/snakelai_out_all/config/slurm-%x_%a.log
#SBATCH --output /shared/projects/elai_most/data/snakelai_out_all/config/slurm-%x_%a.log
#SBATCH --array=[1-22]

module load singularity

## file containing run 1 for each snp set 
wdir=/shared/projects/elai_most/data/snakelai_out_all
dirs=$wdir/local_dosage_dirs.txt

## get one run
l=$SLURM_ARRAY_TASK_ID

elai_dir=$(head -n $l $dirs | tail -1)
echo $elai_dir

## make average results with other runs of the same snp set
singularity exec /shared/projects/vietcaf/rserver/Singularity.R_4-2-0-Rserver_2022.sif Rscript $wdir/config/average_elai.R -e $elai_dir