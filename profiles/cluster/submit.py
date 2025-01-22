#!/usr/bin/env python3
import os
import sys
import logging

from snakemake.utils import read_job_properties

logging.basicConfig(filename='submit_wrapper.log', level=logging.DEBUG)

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

with open(jobscript) as f:
    content = f.read()

logging.debug(content)
# access property defined in the cluster configuration file (Snakemake >=3.6.0)
cluster = job_properties["cluster"]

conda_env = job_properties["cluster"]["conda_env"]

template = f"""#!/usr/bin/env sh

#SBATCH --output=log/%j.out                 # where to store the output ( %j is the JOBID )
#SBATCH --error=log/%j.err                  # where to store error messages

# activate conda environment
. /scicore/home/neher/roemer0001/miniconda3/etc/profile.d/conda.sh
conda activate {conda_env}
export AUGUR_MINIFY_JSON=1
export AUGUR_RECURSION_LIMIT=10000

{content}
"""
logging.debug(template)

with open(jobscript, "w") as f:
    f.write(template)

with open("jobscript.sh", "w") as f:
    f.write(template)

sbatch_command = f"sbatch --time={cluster['time']} --mem={cluster['mem']} --cpus-per-task={cluster['n']} --qos={cluster['qos']} {jobscript}"
logging.debug(sbatch_command)
os.system(sbatch_command)
