mult_ids <- fread('/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/dump/dec_run/merged_ids.rel.id', header=FALSE)

assigned_ids <- fread('/data2/morgante_lab/ukbiobank_projects/GxE_mapping/output/pyse/ethproject.txt')



4120 + 3347 + 746 + 1265 + 15000
idlis <- list.files(path='/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr', pattern='rel.id', full.names = 1)
idlis


id2 <- idlis[grepl('pruned', idlis)]
id2[5] <- gsub('pruned', 'sampled', id2[5])



ids3 <- data.table(path=id2, anc=c('asian', 'black', 'chinese', 'mixed', 'white'))



ids3$pop <- c('SAS', 'AFR','EAS', 'oth', 'EUR')


ids3

hold <- data.table(matrix(c(0,0,0,0), ncol=4))

for(i in 1:nrow(ids3)){
  d1 <- fread(file=unlist(ids3[i,1,with=0]), header=FALSE)
  d1$pop <- ids3[i,3,with=0]
  d1$IID <- paste0(d1$V1, '_',d1$V1)

  hold <- rbind(hold, d1, use.names=FALSE)
}
ids3
i <- 1
ids3[i,1,with=0]


id2[5] <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/sampled_white_ids.rel.id"
fin <- hold[-1,]
names(fin) <- c('FID', 'IID', 'pop', 's')

mids <- assigned_ids[,1:2,]
names(mids) <- c('s', 'rf_anc')



fin
mids


checker <- mids[fin, on=('s')]


checker[,match:=rf_anc==pop]
sum(checker$match)


badmatch <- checker[match==FALSE,]
badmatch
table(badmatch$pop)
table(badmatch$rf_anc)


badmatch
1437 / 24478


unlist(table(badmatch$pop)) / unlist(table(fin$pop)) * 100



otherfix <- badmatch[pop=='oth',]
table(otherfix$rf_anc)



lint <- data.table(names(table(badmatch$pop)),
                   unlist(table(fin$pop)),
                   unlist(table(badmatch$pop)),
                   unlist(table(badmatch$pop)) / unlist(table(fin$pop)) * 100)
outdex <- c(2,4,6)
mint <- lint[,-outdex,with=0]
names(mint) <- c('pop', 'total', 'mismatched', 'percent')
mint$percent <- round(mint$percent, 2)
mint

table(badmatch$rf_anc, by=badmatch$pop)
table(checker$rf_anc, by=checker$pop)




kable(mint, 'simple')
mint2 <- data.table('total', 24478, sum(mint$mismatched), round(sum(mint$mismatched) / 24478 * 100, 2))


mint3 <- rbind(mint, mint2, use.names=FALSE)
mint3




trait   test         P           cutoff   snps
------  -----------  ---------  -------  -----
  DP      SmokingNow   3.16e-06     5.500      2
DP      Poultry      6.31e-06     5.200      3
DP      Cheese       7.24e-05     4.140     31
DP      Tea          1.26e-05     4.900      5
DP      BFP          2.00e-05     4.700      8
SP      TVtime       3.16e-06     5.500      2
SP      Sleepd       1.15e-04     3.940     43
SP      VegCook      2.00e-05     4.700      6
SP      Pork         5.89e-05     4.230     21
SP      Tea          7.00e-05     4.155     25
SP      Coffee       2.00e-05     4.700     10
SP      SleepDev     1.26e-05     4.900      5
PP      TVtime       4.68e-05     4.330     16
PP      Salt         3.16e-05     4.500     11
PP      Tea          3.16e-06     5.500      1
PP      SleepDev     3.16e-04     3.500    114

247+
  125+
  0+
  368+
  77+
  10+
  438

247 / 1265 * 100
125 / 1265 * 100
0 / 1265 * 100
368 / 1265 * 100
77 / 1265 * 100
10 / 1265 * 100
438 / 1265 * 100

20% AFR, 9.8% AMR, 29% EUR, 6% MID, 0.8% other, 34.5% SAS

rule preprocess_ukb_selected_vars:
  vars_csv = "/data2/morgante_lab/data/ukbiobank/variables_processed/UKB_selected_variables.csv",
raw_csv = "/data2/morgante_lab/data/ukbiobank/variables_processed/ukb45105_selected.csv"
output:
  processed_rdata = config["data_dir"] + "/data1_20250501_only_analysis_vars.rds"
params:
  script = config["SCRIPT"] + "/preprocessing_step1_setup_dataset_only_analysis_vars.R"


## filtering dataset ##
rule preprocess_step2:
  input:
  input_rds = config["data_dir"] +"/data1_20250501_only_analysis_vars.rds",
withdrawn = "/data2/morgante_lab/data/ukbiobank/ind_to_remove/withdraw62347_281_20240209.txt",
missing_geno = config["data_dir"] + "/missing_ids.txt"
output:
  filtered_rdata = config["data_dir"] + "/data2_20250501_only_analysis_vars.rds"
params:
  script = config["SCRIPT"] + "/preprocessing_step2_setup_dataset_only_analysis_vars.R"



'/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data1_20250501_only_analysis_vars.rds'



option_list <- list(
  make_option("--input_rds", type = "character"),
  make_option("--withdrawn", type = "character"),
  make_option("--missing_geno", type = "character"),
  make_option("--output_rds", type = "character")
)


opt$input_rds <- '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data1_20250501_only_analysis_vars.rds'
opt$withdrawn <- "/data2/morgante_lab/data/ukbiobank/ind_to_remove/withdraw62347_281_20240209.txt"
opt$missing_geno <- '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/missing_ids.txt'



487979 with genotype
* remove missing DP (31195) or SP (31208)
- 31208 missing DP or SP
456771 remaining
* filter for > 4 SD from mean
- 510 lost
456261



> sum(is.na(dtt$Townsend))
[1] 540
> sum(is.na(dtt$act0_d))
[1] 24241
> sum(is.na(dtt$TVtime))
[1] 25731
> sum(is.na(dtt$sleep_d))
[1] 3423
> sum(is.na(dtt$smoking_now))
[1] 867
> sum(is.na(dtt$veg_cook))
[1] 15522
> sum(is.na(dtt$fish_oily))
[1] 3307
> sum(is.na(dtt$fish_lean))
[1] 3029
> sum(is.na(dtt$meat_proc))
[1] 1671
> sum(is.na(dtt$poultry))
[1] 1524
> sum(is.na(dtt$beef))
[1] 2674
> sum(is.na(dtt$lamb))
[1] 3681
> sum(is.na(dtt$pork))
[1] 3502
> sum(is.na(dtt$cheese))
[1] 12711
> sum(is.na(dtt$salt))
[1] 679
> sum(is.na(dtt$tea))
[1] 15537
> sum(is.na(dtt$alc1))
[1] 1018
> sum(is.na(dtt$waist))
[1] 888
> sum(is.na(dtt$getup))
[1] 1582
> sum(is.na(dtt$coffee))
[1] 34771
> sum(is.na(dtt$smoked_past))
[1] 0
> sum(is.na(dtt$BFP))
[1] 8212


table(dtt$, by=dtt$ethn1_consolidated)
# arth <- rbind(
sum(is.na(dtt$Townsend), by=dtt$ethn1_consolidated),
sum(is.na(dtt$act0_d), by=dtt$ethn1_consolidated),
sum(is.na(dtt$TVtime), by=dtt$ethn1_consolidated),
sum(is.na(dtt$sleep_d), by=dtt$ethn1_consolidated),
sum(is.na(dtt$smoking_now), by=dtt$ethn1_consolidated),
sum(is.na(dtt$veg_cook), by=dtt$ethn1_consolidated),
sum(is.na(dtt$fish_oily), by=dtt$ethn1_consolidated),
sum(is.na(dtt$fish_lean), by=dtt$ethn1_consolidated),
sum(is.na(dtt$meat_proc), by=dtt$ethn1_consolidated),
sum(is.na(dtt$poultry), by=dtt$ethn1_consolidated),
sum(is.na(dtt$beef), by=dtt$ethn1_consolidated),
sum(is.na(dtt$lamb), by=dtt$ethn1_consolidated),
sum(is.na(dtt$pork), by=dtt$ethn1_consolidated),
sum(is.na(dtt$cheese), by=dtt$ethn1_consolidated),
sum(is.na(dtt$salt), by=dtt$ethn1_consolidated),
sum(is.na(dtt$tea), by=dtt$ethn1_consolidated),
sum(is.na(dtt$alc1), by=dtt$ethn1_consolidated),
sum(is.na(dtt$waist), by=dtt$ethn1_consolidated),
sum(is.na(dtt$getup), by=dtt$ethn1_consolidated),
sum(is.na(dtt$coffee), by=dtt$ethn1_consolidated),
sum(is.na(dtt$smoked_past), by=dtt$ethn1_consolidated),
sum(is.na(dtt$BFP), by=dtt$ethn1_consolidated)

arth
library(data.table)

micol <- c('Townsend',
           'act0_d',
           'TVtime',
           'sleep_d',
           'smoking_now',
           'veg_cook',
           'fish_oily',
           'fish_lean',
           'meat_proc',
           'poultry',
           'beef',
           'lamb',
           'pork',
           'cheese',
           'salt',
           'tea',
           'alc1',
           'waist',
           'getup',
           'coffee',
           'smoked_past',
           'BFP')

demo <- setDT(dtt)[, lapply(.SD, function(x) sum(is.na(x))), by = .(ethn1_consolidated), .SDcols=micol]
# G[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = 1:2]


demo <- demo[order(ethn1_consolidated),]
demo
revamp <- data.table(table(dtt$ethn1_consolidated, useNA='always'))
names(revamp) <- c('ethn1', 'SAMPLES')

demo2 <- demo[,-1,with=0]

hold <- data.table(matrix(rep(0,ncol(demo2)), ncol=ncol(demo2)))
names(hold) <- names(demo2)
for (i in 1:nrow(demo2)){
  i <- 1
  newrow <- round(unlist(demo2[i,1:22,with=0]) / unlist(revamp[i,2,with=0]) * 100, digits=1)
  hold <- rbind(hold, newrow)
}

# Asian   Black Chinese   Mixed   White    <NA>
#  9147    7460    1423    2735  429130    6366
