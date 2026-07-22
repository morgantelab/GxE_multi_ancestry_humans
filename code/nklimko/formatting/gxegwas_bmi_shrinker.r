library(data.table)

trait <- 'SP'
env <- 'sleep_dev'
env <- 'cheese'


gxe_shrinker <- function(mytrait, env){
  trait <- toupper(mytrait)

  fullpath <- paste0('/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/',
                     'output/gpcxe/gxegwas_bmi/',
                     trait,
                     '/',
                     env,
                     '.',
                     trait,
                     '0s.glm.linear')

  dataset <- fread(fullpath)

  mitest <- paste0('ADDx',env)

  miset <- dataset[TEST==mitest,]

  micols <- c("#CHROM", 'POS', 'ID', 'BETA', 'SE', 'P')

  fincol <- miset[,micols,with=0]
  fincol$neglogP <- -log10(as.numeric(fincol$P))

  names(fincol)[1:2] <- c('CHR', 'BP')



  outpath <- paste0('/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/gpcxe/tables_gxegwas_bmi/',
                    trait,
                    '/',
                    env,
                    '.rds')

  # print(paste0('outputting to ', outpath))

  saveRDS(fincol, file=outpath)
}

envlist <- c('Townsend', 'act0_d', 'TVtime', 'sleep_d','smoking_now', 'veg_cook', 'fish_oily', 'fish_lean', 'meat_proc', 'poultry', 'beef', 'lamb', 'pork', 'cheese', 'salt', 'tea', 'alc1', 'waist', 'getup', 'coffee', 'smoked_past', 'BFP', 'sleep_dev')

traitlist <- c('dp', 'sp', 'pp')


traitlist <- 'pp'
envlist <- 'Townsend'

for(t in traitlist){
  for(e in envlist){
    gxe_shrinker(mytrait=t, env=e)
  }
}

