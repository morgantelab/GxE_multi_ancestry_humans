#!/bin/bash

#SBATCH --job-name=GE_pcs_PP
#SBATCH --cpus-per-task=15
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=90:30:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/GE_PP_pcrelate_pcs_plink_higher_iter.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/GE_PP_pcrelate_pcs_plink_higher_iter.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

### Set environment variables for threading
source /opt/intel/oneapi/mkl/2023.2.0/env/vars.sh intel64
module load R/4.2.3

export MKL_NUM_THREADS=15

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/creating_vcem_varabs_GE_PP_pcrelate_pcs_plink.R

module unload R/4.2.3
