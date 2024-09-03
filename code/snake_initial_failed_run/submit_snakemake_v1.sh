#!/bin/bash
#
#SBATCH --job-name=snake_genotype_pipeline
#SBATCH --ntasks=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,compute
#SBATCH --time=48:00:00
#SBATCH --mem=6gb
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/snakehead/%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/snakehead/%j.err

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
-n \
--configfile config.yaml \
--latency-wait 30 \
--use-conda \
--nolock \
--rerun-incomplete \
--jobs 100

# Unload R module and deactivate conda environment after completion
#module unload R/4.1.2
conda deactivate

# Log completion message
echo "Snakemake pipeline completed."

