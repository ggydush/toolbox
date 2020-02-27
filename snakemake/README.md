# Snakemake Reference

### Working with Snakemake
1. Download [Snakemake](https://snakemake.readthedocs.io/en/stable/) via `pip install snakemake`
2. If using `qsub`, use `create_snakemake.py` to create `run_snakemake.sh` command to handle job submissions
    * This script uses a qsub wrapper to allow a default configuration to be set for all rules.
    * The jobscript file is used to set up enviornments within qsub
