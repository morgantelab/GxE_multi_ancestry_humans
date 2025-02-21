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