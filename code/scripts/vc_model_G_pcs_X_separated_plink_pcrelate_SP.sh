#!/bin/bash

#SBATCH --job-name=sp_plink_pcrelate_X_sep
#SBATCH --cpus-per-task=5
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=90:30:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/X_separated_model_pcs_plink_sp_pcrelate_25.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/X_separated_model_pcs_plink_sp_pcrelate_25.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

###Load module
source /opt/intel/oneapi/mkl/2023.2.0/env/vars.sh intel64
module load R/4.2.3

export MKL_NUM_THREADS=5
Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/vc_model_G_pcs_X_separated_plink_pcrelate_SP.R

module unload R/4.2.3
