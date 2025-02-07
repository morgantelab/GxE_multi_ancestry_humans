#!/bin/bash
#
#SBATCH --job-name=gwas_dp
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute,fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=1:00:00
#SBATCH --mem=100G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/test_gwas_plink_act0_d_train.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm/test_gwas_plink_act0_d_train.err
#SBATCH --mail-type=all
#SBATCH --mail-user=kgoda@clemson.edu

##
/data2/morgante_lab/kgoda/software/plink2 --bfile /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv \
--keep /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/Train_Set_1.txt \
--pheno /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/covar_pheno.txt --pheno-name DP0s \
--covar /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/covar_pheno.txt --covar-name age,age2,sex,X1-X10,act0_d \
--glm interaction \
--parameters 1-15,29 \
--out /data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/plink_gwas_snp_act0_d_train
