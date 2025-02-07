#!/bin/bash

#SBATCH --job-name=eigen_GE
#SBATCH --cpus-per-task=10
#SBATCH --partition=compute
#SBATCH --time=50:30:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/eigen_GE_25.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/eigen_GE_25.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

module load R/4.1.2

Rscript /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/creating_eigen_GE_hadamard_prods.R

module unload R/4.1.2
