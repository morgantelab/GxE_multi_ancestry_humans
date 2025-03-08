#!/bin/bash
#
#SBATCH --job-name=snake_genotype_pipeline
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --time=10-00:00:00
#SBATCH --mem=1gb
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/log/Emat_S_A_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/log/Emat_S_A_%j.err

# Log message to confirm script starts running
echo "Starting Snakemake pipeline..."

# Create necessary directories for logs
#mkdir -p /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/snakehead

# Navigate to the working directory where the Snakefile and config.yaml are located
cd /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code

# Source conda/mamba and activate the snakemake environment
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh
conda activate snakemake

# Load the R module
#module load R/4.1.2

# Set environment variables to limit threading in linear algebra libraries
#export OPENBLAS_NUM_THREADS=1
#export OMP_NUM_THREADS=1

# Include personal R packages
#export R_LIBS=/data2/morgante_lab/kgoda/software/R/x86_64-pc-linux-gnu-library/4.1

# Execute Snakemake
snakemake \
-s snakefile \
--configfile config.yaml \
--latency-wait 120 \
--use-conda \
--keep-going \
--rerun-incomplete \
--profile slurm

# Test/Dry version
#snakemake \
#-p \
#-n \
#-s snakefile.yaml \
#--configfile config.yaml \
#--profile slurm

#--dag | display | dot
#-p -n \

# Unload R module and deactivate conda environment after completion
#module unload R/4.1.2
conda deactivate

# Log completion message
echo "Snakemake pipeline completed."

