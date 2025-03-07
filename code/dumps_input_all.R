## rule alls of snakefile ##

rule all:
  input:
  expand(config["gwas_dir"] + "/summary_hits_{pval}.txt", pval=config["PVAL"])

rule all:
  input:
  expand(config["gwas_dir"] + "/Train_Set_{foldnum}.txt", foldnum=config["FOLDNUM"]),
  expand(config["gwas_dir"] + "/Test_Set_{foldnum}.txt", foldnum=config["FOLDNUM"]),
  expand(config["gwas_dir"] + "/Train_exclude_{ancestry}.txt", ancestry=config["Ancestry"]),
  expand(config["gwas_dir"] + "/Test_{ancestry}.txt", ancestry=config["Ancestry"])

rule all:
  input:
  expand(config["gwas_dir"] + "/Hadamard_GRM_Fold{foldnum}_{bp}_{pval}.RData", foldnum=config["FOLDNUM"], bp=config["BP"], pval=config["PVAL"]),
expand(config["gwas_dir"] + "/Hadamard_GRM_Ancestry{ancestry}_{bp}_{pval}.RData", ancestry=config["Ancestry"], bp=config["BP"], pval=config["PVAL"])

rule all:
  input:
  expand(config["sample_ids_dir"] + "/" + "Emat_with_sex_age.RDS"),
  expand(config["sample_ids_dir"] + "/" + "Emat_eigen_with_sex_age.RDS")

rule all:
  input:
  expand(config["snp_list_dir"] + "/" + "{ancestry}_filtered_all_chr{chr}.snplist", chr=config["CHR"], ancestry=config["Ancestry"])

rule all:
  input:
  expand(config["snp_list_dir"] + "/" + "common_snps_all_ancestries.txt")

rule all:
  input:
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_G.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_E.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_G_E.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_G_E_GE.csv", bp=config["BP"], ancestry=config["Ancestry"])

rule all:
  input:
  expand(config["model_dir"] + "/" + "VCEm_{bp}_{grm}_G.csv", bp=config["BP"], grm=config["GRM"]),
  expand(config["model_dir"] + "/" + "VCEm_{bp}_{grm}_G_pcs.csv", bp=config["BP"], grm=config["GRM"]),
  expand(config["model_dir"] + "/" + "VCEm_{bp}_{grm}_G_pcs_X_separated.csv", bp=config["BP"], grm=config["GRM"]),
  expand(config["model_dir"] + "/" + "VCEm_{bp}_{grm}_G_E.csv", bp=config["BP"], grm=config["GRM"]),
  #expand(config["model_dir"] + "/" + "VCEm_{bp}_{grm}_run_G_E_pcrelate_pcs_plink.csv", bp=config["BP"], grm=config["GRM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_Fold_{foldnum}_X1.csv", bp=config["BP"], foldnum=config["FOLDNUM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_Fold_{foldnum}_X1_X2.csv", bp=config["BP"], foldnum=config["FOLDNUM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_Fold_{foldnum}_X1_X2_G.csv", bp=config["BP"], foldnum=config["FOLDNUM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2.csv",bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_G.csv",bp=config["BP"], ancestry=config["Ancestry"])

rule all:
  input:
  expand(config["model_dir"] + "/" + "VCEm_{bp}_{grm}_run_G_E_pcrelate_pcs_plink.csv", bp=config["BP"], grm=config["GRM"]),
  expand(config["model_dir"] + "/" + "varabs_{bp}_{grm}_GE.csv", bp=config["BP"], grm=config["GRM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_Fold_{foldnum}_X1_X2_E.csv", bp=config["BP"], foldnum=config["FOLDNUM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_Fold_{foldnum}_X1_X2_G_E.csv", bp=config["BP"], foldnum=config["FOLDNUM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_Fold_{foldnum}_X1_X2_G_E_GE.csv", bp=config["BP"], foldnum=config["FOLDNUM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_E.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_G_E.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_G_E_GE.csv", bp=config["BP"], ancestry=config["Ancestry"])

rule all:
  input:
  expand(config["model_dir"] + "/" + "varabs_{bp}_{grm}_GE.csv", bp=config["BP"], grm=config["GRM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_Fold_{foldnum}_X1_X2_E.csv", bp=config["BP"], foldnum=config["FOLDNUM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_Fold_{foldnum}_X1_X2_G_E.csv", bp=config["BP"], foldnum=config["FOLDNUM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_Fold_{foldnum}_X1_X2_G_E_GE.csv", bp=config["BP"], foldnum=config["FOLDNUM"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_E.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_G_E.csv", bp=config["BP"], ancestry=config["Ancestry"]),
  expand(config["model_dir"] + "/" + "PREDs_{bp}_ethn_{ancestry}_X1_X2_G_E_GE.csv", bp=config["BP"], ancestry=config["Ancestry"])





