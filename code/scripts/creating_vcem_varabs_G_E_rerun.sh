#!/bin/bash

#SBATCH --job-name=rerun_G_E
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=90:30:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/rerun_G_E.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/re_run_G_E.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

### Set environment variables for threading
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

module load R/4.1.2

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/creating_vcem_varabs_G_E.R

module unload R/4.1.2