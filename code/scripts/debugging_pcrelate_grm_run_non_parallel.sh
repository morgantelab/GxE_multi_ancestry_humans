#!/bin/bash

#SBATCH --job-name=pcrelate_test
#SBATCH --cpus-per-task=55
#SBATCH --partition=compute,fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=50:30:00
#SBATCH --mem=900G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/test_pcrelate_single.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/test_pcrelate_single.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

source /opt/intel/oneapi/mkl/2023.2.0/env/vars.sh intel64
module load R/4.2.3

export MKL_NUM_THREADS=55

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/debugging_pcrelate_grm_run_non_parallel.R

module unload R/4.2.3