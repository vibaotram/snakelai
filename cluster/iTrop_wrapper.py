#!/usr/bin/env python3

'''
wrap sbatch jobs based on Snakemake config file
'''
import os
import sys
from snakemake.utils import read_job_properties
import re
# import subprocess

jobscript = sys.argv[-1]

job_properties = read_job_properties(jobscript)


jobname = job_properties['params']['logname']
job_name = '--job-name ' + jobname

# jobid = job_properties['jobid']
threads = job_properties['threads']
cpus_per_task = '--cpus-per-task ' + str(threads)

ntasks = '--ntasks 1'


logdir = os.path.join(job_properties['params']['logdir'], "slurm")
os.makedirs(logdir, exist_ok=True)
try:
    log = job_properties['params']['logname']
    # log = os.path.splitext(log)[0]
except IndexError:
    log = rule
output = f'--output {logdir}/{log}_%j.log'
error = f'--error {logdir}/{log}_%j.log'

try:
    mem_gb = job_properties['resources']['mem_gb']
    mem = f'--mem {mem_gb}G'
except (IndexError, KeyError):
    mem = "--mem 10G"

sbatch_params = ' '.join(['sbatch --parsable -p long', job_name, cpus_per_task, mem, ntasks, output, error])

dep_jobid = sys.argv[1:-1]
if not any(re.match('\d+', j) for j in dep_jobid): # "normal"-submit
    dependencies = ''
else: # if immediate-submit, for slurm only
    dependencies = ' --dependency=afterok:' + ':'.join(dep_jobid)

sbatch_params += dependencies

cmdline = " ".join([sbatch_params, jobscript])

# jobscript.replace("\n", "\necho -e\"sbatch parameters:\n\"{}\"\"".format(sbatch), 1)
with open(jobscript, "r") as j:
    scripts = j.readlines()
scripts.insert(1, "echo -e \"# Submit command-line: \"{}\"\"\n".format(cmdline))
scripts.insert(2, "echo -e \"# Job running on node: $SLURM_JOB_NODELIST\"\n")
scripts.insert(3, "echo -e \"\n\"")
with open(jobscript, "w") as j:
    j.writelines(scripts)

os.system(cmdline)
# print(cmdline)
# sbatch --job-name {cluster.job-name} --partition {cluster.partition} --account {cluster.account} --cpus-per-task {cluster.cpus-per-task} --output {cluster.output} --error {cluster.error} snakejob.sh
