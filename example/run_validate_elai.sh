#!/bin/bash
#SBATCH --job-name benchmark_elai
#SBATCH --partition long
#SBATCH --cpus-per-task 11
#SBATCH --output /shared/projects/elai_most/benchmark_elai/config/slurm-%x.log
#SBATCH --error /shared/projects/elai_most/benchmark_elai/config/slurm-%x.log

module load singularity
module load conda
module load snakemake/6.5.0

cd /shared/projects/vietcaf/scripts/snakelai
source venv/bin/activate

snakemake --nolock --use-conda --use-singularity \
--singularity-args " -B /shared/home/baotram,/shared/projects/vietcaf,/shared/projects/elai_most/benchmark_elai" \
--cores --jobs -p --verbose \
--latency-wait 60 --keep-going --rerun-incomplete \
--cluster "python3 cluster/iTrop_wrapper.py" \
--cluster-status "python3 cluster/iTrop_status.py" \
--configfile "/shared/projects/elai_most/benchmark_elai/config/config.yaml"
