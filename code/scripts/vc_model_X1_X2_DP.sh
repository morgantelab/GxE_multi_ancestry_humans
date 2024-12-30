#!/bin/bash
#SBATCH --job-name=Xs_DP
#SBATCH --cpus-per-task=10
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=30:30:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/X1_X2_DP.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/X1_X2_DP.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

###Load module
source /opt/intel/oneapi/mkl/2023.2.0/env/vars.sh intel64
module load R/4.2.3

export MKL_NUM_THREADS=10

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/vc_model_just_X1_X2_DP.R

module unload R/4.2.3