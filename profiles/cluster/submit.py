#!python

#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# access property defined in the cluster configuration file (Snakemake >=3.6.0)
cluster = job_properties["cluster"]

conda_env = job_properties["cluster"]["conda_env"]

with open(jobscript,'a') as f:
    f.write(f"conda activate {conda_env};")
    f.write("{exec_job}")

os.system(f"sbatch --time={cluster['time']} --mem={cluster['mem']} --cpus-per-task={cluster['n']} --qos={cluster['qos']} {jobscript}")