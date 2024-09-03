#!/bin/bash
mkdir -p ../output/{log,logs_slurm} | sbatch ./submit_snakemake.sh
