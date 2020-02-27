#!/usr/bin/env python

import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

defaults = {"mem": 8, "runtime": 2}
params = {**defaults, **job_properties["resources"]}

qsub_cmd = (
    f'qsub -l h_vmem={params["mem"]}G '
    f'-l h_rt={params["runtime"]}:00:00 '
    f'-o logs/{job_properties["rule"]}/ '
    f'-N {job_properties["rule"]} '
    f"-cwd -j y -sync y "
    f"{jobscript}"
)

qsub_cmd += ' | tail -2 | cut -d " " -f 3'
os.system(qsub_cmd)
