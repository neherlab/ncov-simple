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

print(jobscript)
print("hello test")

command = jobscript.splitlines()[-1]

template = f"""
#!/bin/sh

#SBATCH --output=log/%j.out                 # where to store the output ( %j is the JOBID )
#SBATCH --error=log/%j.err                  # where to store error messages

# activate conda environment
source /scicore/home/neher/roemer0001/miniconda3/etc/profile.d/conda.sh
conda activate {conda_env}
export AUGUR_MINIFY_JSON=1
export AUGUR_RECURSION_LIMIT=10000

{command}
"""

os.system(f"sbatch --time={cluster['time']} --mem={cluster['mem']} --cpus-per-task={cluster['n']} --qos={cluster['qos']} {template}")