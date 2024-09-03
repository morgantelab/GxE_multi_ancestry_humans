#!/bin/bash
mkdir -p ./{../output,../output/logs_slurm} | sbatch ./submit_snakemake.sh
