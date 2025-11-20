# GxE Multi-Ancestry Analysis

This project performs genome-wide gene-by-environment (GxE) interaction analysis of blood pressure traits using UK Biobank data across multiple ancestries.

## Overview

- All analysis steps are orchestrated using a `Snakefile`, which defines the full computational workflow including data preprocessing, modeling, GxE-GWAS, and post-analysis.
- The `Snakefile` contains the complete set of rules and dependencies to ensure reproducibility of all analytical steps.
- The 'snakefile' can be found in /data2/morgante_lab/ukbiobak_projects/GxE_multi_ancestry/code/snakefile.
- the corresponding config and initiator files are also in same dir. slurm folder has the snakemake config file.

## Scripts

All R and shell scripts used in the `Snakefile` are stored in the following directory:

/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/code/scripts/snakemake_scripts/ 

Each script corresponds to a specific analysis rule and is modularized for reusability and clarity.

Scripts for plots used in Manuscript are in scripts/plot_scripts/manuscript_plots/ folder.

All additional scripts of analysis ran but not used, plots generated but not used in manuscript can be found in code/scripts/processing_scripts/

## Usage

To run the analysis:

"bash initiator.sh"

## Runs logs

All log files for each snakefile run can be found in output/log folder.
All rule ran log files can be found in output/logs_slurm folder.

## Data and Results

All data generated from the snakefile analysis runs can be found in /data folder.

# Paper

This repository contains code and data resources to accompany our research paper:

> Goda, K., Klimkowski Arango, N., Tiezzi, F., Mackay, T. F. C., & Morgante, F. (2025).
> Gene-environment interactions contribute to blood pressure variation across global populations.
> *medRxiv* 2025.07.02.25330727.
> https://doi.org/10.1101/2025.07.02.25330727
