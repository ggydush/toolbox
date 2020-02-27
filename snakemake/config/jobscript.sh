#!/bin/bash
# properties = {properties}
. /broad/software/free/Linux/redhat_6_x86_64/pkgs/anaconda3_5.0.1/etc/profile.d/conda.sh;
conda deactivate && conda activate /xchip/bloodbiopsy/envs/bloodbiopsy;
export PYTHONPATH=$PYTHONPATH:/xchip/bloodbiopsy/rhoades/scripts/bbutils;
source /broad/software/scripts/useuse;
reuse -q UGER;
reuse -q R-3.4;
reuse -q Bowtie2;
{exec_job}
