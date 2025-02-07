#!/bin/bash

#SBATCH --job-name=ks_test
#SBATCH --cpus-per-task=10
#SBATCH --partition=compute,fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=3:30:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/ks_test.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/ks_test.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

source /opt/intel/oneapi/mkl/2023.2.0/env/vars.sh intel64
module load R/4.2.3

export MKL_NUM_THREADS=1

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/ks_test_pvalue_uniformity.R

module unload R/4.2.3