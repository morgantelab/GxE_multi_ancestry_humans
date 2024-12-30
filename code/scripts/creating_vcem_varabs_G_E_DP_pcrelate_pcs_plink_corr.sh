#!/bin/bash

#SBATCH --job-name=G_E_pcs_DP
#SBATCH --cpus-per-task=10
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=90:30:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/G_E_DP_pcrelate_pcs_plink_2.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/G_E_DP_pcrelate_pcs_plink_2.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

### Set environment variables for threading
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

module load R/4.1.2

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/creating_vcem_varabs_G_E_DP_pcrelate_pcs_plink_corr.R

module unload R/4.1.2