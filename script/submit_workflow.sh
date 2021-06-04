#!/bin/bash
#SBATCH --job-name elai_TR6
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log

module load system/Miniconda3/1.0

source activate snakemake

snakemake --nolock --use-singularity --singularity-args "--nv " --use-conda --cores -p --verbose \
--latency-wait 60 --keep-going --restart-times 2 --rerun-incomplete \
--cluster "python3 cluster/iTrop_wrapper.py" \
--cluster-status "python3 cluster/iTrop_status.py" \
--configfile "/data3/projects/vietcaf/baotram/scripts/robusta_vn/test_10_clones/elai_tr6/configs/config_elai_tr6_even_snps.yaml"
