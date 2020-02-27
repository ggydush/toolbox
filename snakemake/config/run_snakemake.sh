snakemake --cluster-sync $PWD/qsub_wrapper.py --jobscript $PWD/jobscript.sh -j 50 --latency-wait 120 --restart-times 3
