#!/bin/bash

#SBATCH --job-name=pred_asian
#SBATCH --cpus-per-task=10
#SBATCH --partition=compute,fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=24:00:00
#SBATCH --mem=150G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/pred_asian_GEselect.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/pred_asian_GEselect.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

source /opt/intel/oneapi/mkl/2023.2.0/env/vars.sh intel64
module load R/4.2.3

export MKL_NUM_THREADS=10

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/pred_model_ethn_X1_X2_G_E_GEselect_1e3_check.R

module unload R/4.2.3
