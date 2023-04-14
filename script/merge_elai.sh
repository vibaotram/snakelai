#!/bin/bash
#SBATCH --job-name merge_elai
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --error /shared/projects/elai_most/data/snakelai_out_all/config/slurm-%x_%a.log
#SBATCH --output /shared/projects/elai_most/data/snakelai_out_all/config/slurm-%x_%a.log
#SBATCH --array=[8]

module load singularity

## chrom number
range=$(printf %02d $SLURM_ARRAY_TASK_ID)

## merge elai results for the chrom
singularity exec /shared/projects/vietcaf/rserver/Singularity.R_4-2-0-Rserver_2022.sif Rscript /shared/projects/elai_most/data/snakelai_out_all/config/merge_elai.R $range
