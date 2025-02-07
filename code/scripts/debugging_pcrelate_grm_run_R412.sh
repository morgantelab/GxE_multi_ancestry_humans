#!/bin/bash

#SBATCH --job-name=pcrelate_test_412
#SBATCH --cpus-per-task=55
#SBATCH --partition=compute,fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=50:30:00
#SBATCH --mem=900G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/test_pcrelate_R_412.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/test_pcrelate_R_412.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

module load R/4.1.2

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/debugging_pcrelate_grm_run_R412.R

module unload R/4.1.2
