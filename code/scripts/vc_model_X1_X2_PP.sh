#!/bin/bash
#SBATCH --job-name=Xs_PP
#SBATCH --cpus-per-task=5
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=30:30:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/X1_X2_PP_25.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/X1_X2_PP_25.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

###Load module
source /opt/intel/oneapi/mkl/2023.2.0/env/vars.sh intel64
module load R/4.2.3

export MKL_NUM_THREADS=5

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/vc_model_X1_X2_PP.R

module unload R/4.2.3
