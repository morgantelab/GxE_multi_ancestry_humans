#!/bin/bash
#SBATCH --job-name=test_R_10
#SBATCH --cpus-per-task=10
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=90:30:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/test_R_cpus10.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/test_R_cpus10.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

###Load module
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

export MKL_NUM_THREADS=10

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/test_R_versions_cpus1.R

module unload R/4.2.3