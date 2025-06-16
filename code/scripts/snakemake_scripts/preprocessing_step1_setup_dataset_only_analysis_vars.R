set.seed(1123)

rm(list=ls()); gc()

library(optparse)
library(data.table)

option_list <- list(
  make_option(c("--vars_csv"), type = "character", help = "CSV with variable metadata"),
  make_option(c("--raw_csv"), type = "character", help = "Raw UKB data file"),
  make_option(c("--out_rdata"), type = "character", help = "Output RDS file path")
)
opt <- parse_args(OptionParser(option_list = option_list))

vars_csv <- opt$vars_csv
raw_csv <- opt$raw_csv
out_rdata <- opt$out_rdata

# Load variable description
dat_codes <- fread(vars_csv, data.table = FALSE)
dat_codes <- subset(dat_codes, to_be_kept == "Yes")

# Load raw data
dat <- fread(raw_csv, data.table = FALSE)
colnames(dat) <- paste("V", colnames(dat), sep = '')




# 1230	6177-0.0	226928	Categorical (multiple)





#############
##### P #####
#############

dat$DP0 <- dat$DP0a <- dat$'V4079-0.0'
dat$SP0 <- dat$SP0a <- dat$'V4080-0.0'
dat$PP0 <- dat$SP0-dat$DP0

phen <- c("diast_blood_00",
          "sist_blood_00",
          "pulse_pressure_00")
phen <- c("DP0", "SP0", "PP0", "DP0a", "SP0a")
data.frame(phen,phen%in%colnames(dat))







#############
##### I #####
#############

# ID  # eid
dat$ID <- dat$Veid

# Sex # V31-0.0
dat$Sex_SI <- dat$'V31-0.0'

# Year of Birth # V34-0.0
dat$YOB <- dat$'V34-0.0'

# Month of Birth # V52-0.0
dat$MOB <- dat$'V52-0.0'

# Date of phenotyping # V53-0.0
dat$DOF <- dat$'V53-0.0'

# Phenotyping center # V54-0.0
dat$COF <- dat$'V54-0.0'

# Ethnicity based on UK Biobank Data-Field 21000
dat$ethn1 <- dat$'V21000-0.0'

# Replace 'Do not know' (-1) and 'Prefer not to answer' (-3) with NA
dat$ethn1 <- replace(dat$ethn1, dat$ethn1 %in% c(-1, -3), NA)

# Grouping ethnicity categories
dat$ethn1_whbri <- ifelse(dat$ethn1 == 1001, 1, 0)  # White British
dat$ethn1_white <- ifelse(dat$ethn1 %in% c(1, 1001, 1002, 1003), 1, 0)  # Any White background
dat$ethn1_mixed <- ifelse(dat$ethn1 %in% c(2, 2001, 2002, 2003, 2004), 1, 0)  # Mixed background
dat$ethn1_asian <- ifelse(dat$ethn1 %in% c(3, 3001, 3002, 3003, 3004), 1, 0)  # Asian or Asian British
dat$ethn1_black <- ifelse(dat$ethn1 %in% c(4, 4001, 4002, 4003), 1, 0)  # Black or Black British
dat$ethn1_chinese <- ifelse(dat$ethn1 == 5, 1, 0)  # Chinese
dat$ethn1_other <- ifelse(dat$ethn1 == 6, 1, 0)  # Other ethnic group


# Age phenotyping, recruitment # V21003-0.0 V21022-0.0
dat$AOP <- dat$'V21003-0.0'
dat$AOR <- dat$'V21022-0.0'

Ivars <- c("ID", "Sex_SI", "YOB", "MOB", "DOF", "COF", "ethn1", "ethn1_white", "ethn1_whbri", "ethn1_mixed", "ethn1_asian", "ethn1_black", "ethn1_chinese", "ethn1_other", "AOP", "AOR")
data.frame(Ivars,Ivars%in%colnames(dat))

#############
##### E #####
#############


##### Body #####

# Waist circumference # V48-0.0
dat$waist <- dat$'V48-0.0'

# Body fat percentage # V23099-0.0
dat$BFP <- dat$'V23099-0.0'

# Basal metabolic rate # V23105-0.0
dat$BMR <- dat$'V23105-0.0'

Ebody <- c("waist", "BFP", "BMR")
data.frame(Ebody,Ebody%in%colnames(dat))

##### Illness #####

# Cancer year/age first occurred # table(dat$'V87-0.0')

# Number of cancer, non-cancer, medications self-reported # V134-0.0 V135-0.0 V137-0.0
dat$cancer1 <- ifelse(dat$'V134-0.0'==0, 0, 1)
dat$cancer2 <- dat$'V2453-0.0'; dat$cancer2 <- replace(dat$cancer2, dat$'V2453-0.0'<0, NA)
dat$disease <- ifelse(dat$'V135-0.0'==0, 0, 1)
dat$medicat1 <- dat$'V137-0.0'
dat$medicat2 <- dat$'V2492-0.0'; dat$medicat2 <- replace(dat$medicat2, dat$'V2492-0.0'<0, NA)

# medications for pain, constipation, heartburn #V6154-0.0
dat$medicat3 <- 0
dat$medicat3 <- replace(dat$medicat3, dat$'V6154-0.0'%in%c(-7), 0)
dat$medicat3 <- replace(dat$medicat3, dat$'V6154-0.0'%in%c(-1, -3), NA)
dat$medicat3 <- replace(dat$medicat3, dat$'V6154-0.0'%in%c(1:6), 1)

# Falls in the last year # V2296-0.0
dat$falls <- 0
dat$falls <- replace(dat$falls, dat$'V2296-0.0'%in%c(-3), NA)
dat$falls <- replace(dat$falls, dat$'V2296-0.0'%in%c(1), 0)
dat$falls <- replace(dat$falls, dat$'V2296-0.0'%in%c(2,3), 1)



# Medication for cholesterol, blood pressure or diabetes

# 6153-0.0	270845	Categorical (multiple)	Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones
dat$mediblood1 <- 0
dat$mediblood1 <- replace(dat$mediblood1, dat$'V6153-0.0'%in%c(-7,1,3,4,5), -1)
dat$mediblood1 <- replace(dat$mediblood1, dat$'V6153-0.0'%in%c(-1, -3), NA)
dat$mediblood1 <- replace(dat$mediblood1, dat$'V6153-0.0'%in%c(2), 1)

# 6153-0.1	270845	Categorical (multiple)	Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones
dat$mediblood3 <- 0
dat$mediblood3 <- replace(dat$mediblood3, dat$'V6153-0.1'%in%c(-7,1,3,4,5), -1)
dat$mediblood3 <- replace(dat$mediblood3, dat$'V6153-0.1'%in%c(-1, -3), NA)
dat$mediblood3 <- replace(dat$mediblood3, dat$'V6153-0.1'%in%c(2), 1)

# 6153-0.2	270845	Categorical (multiple)	Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones
dat$mediblood5 <- 0
dat$mediblood5 <- replace(dat$mediblood5, dat$'V6153-0.2'%in%c(-7,1,3,4,5), -1)
dat$mediblood5 <- replace(dat$mediblood5, dat$'V6153-0.2'%in%c(-1, -3), NA)
dat$mediblood5 <- replace(dat$mediblood5, dat$'V6153-0.2'%in%c(2), 1)

# 6153-0.3	270845	Categorical (multiple)	Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones
dat$mediblood6 <- 0
dat$mediblood6 <- replace(dat$mediblood6, dat$'V6153-0.3'%in%c(-7,1,3,4,5), -1)
dat$mediblood6 <- replace(dat$mediblood6, dat$'V6153-0.3'%in%c(-1, -3), NA)
dat$mediblood6 <- replace(dat$mediblood6, dat$'V6153-0.3'%in%c(2), 1)

# 6177-0.0	226928	Categorical (multiple)
dat$mediblood2 <- 0
dat$mediblood2 <- replace(dat$mediblood2, dat$'V6177-0.0'%in%c(-7,1,3,4,5), -1)
dat$mediblood2 <- replace(dat$mediblood2, dat$'V6177-0.0'%in%c(-1, -3), NA)
dat$mediblood2 <- replace(dat$mediblood2, dat$'V6177-0.0'%in%c(2), 1)

# 6177-0.1	226928	Categorical (multiple)
dat$mediblood4 <- 0
dat$mediblood4 <- replace(dat$mediblood4, dat$'V6177-0.1'%in%c(-7,1,3,4,5), -1)
dat$mediblood4 <- replace(dat$mediblood4, dat$'V6177-0.1'%in%c(-1, -3), NA)
dat$mediblood4 <- replace(dat$mediblood4, dat$'V6177-0.1'%in%c(2), 1)

# 6177-0.2	226928	Categorical (multiple)
dat$mediblood7 <- 0
dat$mediblood7 <- replace(dat$mediblood7, dat$'V6177-0.2'%in%c(-7,1,3,4,5), -1)
dat$mediblood7 <- replace(dat$mediblood7, dat$'V6177-0.2'%in%c(-1, -3), NA)
dat$mediblood7 <- replace(dat$mediblood7, dat$'V6177-0.2'%in%c(2), 1)


# diabetes diagnoses # V2443-0.0
dat$diab <- dat$'V2443-0.0'; dat$diab <- replace(dat$diab, dat$diab<0, NA)

# vascular problems #
dat$vascul <- 0
dat$vascul <- replace(dat$vascul, dat$'V6150-0.0'%in%c(-3), NA)
dat$vascul <- replace(dat$vascul, dat$'V6150-0.0'%in%c(-7), 0)
dat$vascul <- replace(dat$vascul, dat$'V6150-0.0'%in%c(1,2,3,4), 1)

# doctor diagnosis # V6152-0.0
dat$doctor <- 0
dat$doctor <- replace(dat$doctor, dat$'V6152-0.0'%in%c(-3), NA)
dat$doctor <- replace(dat$doctor, dat$'V6152-0.0'%in%c(-7), 0)
dat$doctor <- replace(dat$doctor, dat$'V6152-0.0'%in%c(1:9), 1)

# Non cancer illness, treatments, age at #V20002-0.0 V20003-0.0 V20009-0.0
dat$illness <- dat$'V20002-0.0'
dat$treat <- dat$'V20003-0.0'
dat$ageill <- dat$'V20009-0.0'
dat$ageill <- replace(dat$ageill, dat$ageill%in%c(-1,-3), NA)

# FEV1, FVC # V20150-0.0 V20151-0.0
dat$FEV1 <- dat$'V20150-0.0'
dat$FVC <- dat$'V20151-0.0'

Eill <- c("cancer1", "cancer2", "disease", "medicat1", "medicat2", "medicat3", "mediblood1", "mediblood2", "mediblood3", "mediblood4", "mediblood5", "mediblood6", "mediblood7", "falls", "diab", "vascul", "doctor", "illness", "treat", "ageill", "FEV1", "FVC")
data.frame(Eill,Eill%in%colnames(dat))







##### Life #####
# Townsend deprivation index at recruitment # V189-0.0
dat$Townsend <- dat$'V189-0.0'

Elife <- c("Townsend")
data.frame(Elife,Elife%in%colnames(dat))

##### Physical activity #####


# Number of days/week walked 10+ minutes
dat$walk_d <- dat$'V864-0.0'
dat$walk_d <- replace(dat$walk_d, dat$walk_d<0, NA)

# Number of days/week of moderate physical activity 10+ minutes
dat$act0_d <- dat$'V884-0.0'
dat$act0_d <- replace(dat$act0_d, dat$act0_d<0, NA)

#Number of days/week of vigorous physical activity 10+ minutes
dat$act1_d <- dat$'V904-0.0'
dat$act1_d <- replace(dat$act1_d, dat$act1_d<0, NA)

# Time spent watching TV
dat$TVtime <- dat$'V1070-0.0'
dat$TVtime <- replace(dat$TVtime, dat$TVtime<0, NA)
dat$TVtime <- replace(dat$TVtime, dat$TVtime==(-10), 0)

# Time spent using computer
dat$PCtime <- dat$'V1080-0.0'
dat$PCtime <- replace(dat$PCtime, dat$PCtime<0, NA)
dat$PCtime <- replace(dat$PCtime, dat$PCtime==(-10), 0)

# Time spent driving
dat$DRtime <- dat$'V1090-0.0'
dat$DRtime <- replace(dat$DRtime, dat$DRtime<0, NA)
dat$DRtime <- replace(dat$DRtime, dat$DRtime==(-10), 0)

Ephys <- c("walk_d", "act0_d", "act1_d", "TVtime", "PCtime", "DRtime")
data.frame(Ephys,Ephys%in%colnames(dat))





##### Resting habit #####

# Sleep duration
dat$sleep_d <- dat$'V1160-0.0'
dat$sleep_d <- replace(dat$sleep_d, dat$sleep_d%in%c(-3,-1), NA)

# Getting up in the morning, 6 classes
dat$getup <- dat$'V1170-0.0'
dat$getup <- replace(dat$getup, dat$getup%in%c(-3,-1), NA)

Esleep <- c("sleep_d", "getup")
data.frame(Esleep,Esleep%in%colnames(dat))









##### Smoking #####

# Tobacco smoking, current
dat$smoking_now <- ifelse(dat$'V1239-0.0'>0, 1, 0)
dat$smoking_now <- replace(dat$smoking_now, dat$'V1239-0.0'<0, NA)

# Tobacco smoking, past
dat$smoked_past <- ifelse(dat$'V1249-0.0'%in%c(1,2,3), 1, 0)

Esmoke <- c("smoking_now", "smoked_past")
data.frame(Esmoke,Esmoke%in%colnames(dat))







##### Food #####

# Cooked vegetable intake
dat$veg_cook <- dat$'V1289-0.0'
dat$veg_cook <- replace(dat$veg_cook, dat$veg_cook<0, NA)
dat$veg_cook <- replace(dat$veg_cook, dat$veg_cook==(-10), 0)

# Oily fish intake
dat$fish_oily <- dat$'V1329-0.0'
dat$fish_oily <- replace(dat$fish_oily, dat$fish_oily<0, NA)
dat$fish_oily <- replace(dat$fish_oily, dat$fish_oily==(-10), 0)

# Non-Oily fish intake
dat$fish_lean <- dat$'V1339-0.0'
dat$fish_lean <- replace(dat$fish_lean, dat$fish_lean<0, NA)
dat$fish_lean <- replace(dat$fish_lean, dat$fish_lean==(-10), 0)

# Processed meat intake
dat$meat_proc <- dat$'V1349-0.0'
dat$meat_proc <- replace(dat$meat_proc, dat$meat_proc<0, NA)
dat$meat_proc <- replace(dat$meat_proc, dat$meat_proc==(-10), 0)

# Poultry intake
dat$poultry <- dat$'V1359-0.0'
dat$poultry <- replace(dat$poultry, dat$poultry<0, NA)
dat$poultry <- replace(dat$poultry, dat$poultry==(-10), 0)

# Beef intake
dat$beef <- dat$'V1369-0.0'
dat$beef <- replace(dat$beef, dat$beef<0, NA)
dat$beef <- replace(dat$beef, dat$beef==(-10), 0)

# Lamb intake
dat$lamb <- dat$'V1379-0.0'
dat$lamb <- replace(dat$lamb, dat$lamb<0, NA)
dat$lamb <- replace(dat$lamb, dat$lamb==(-10), 0)

# Pork intake
dat$pork <- dat$'V1389-0.0'
dat$pork <- replace(dat$pork, dat$pork<0, NA)
dat$pork <- replace(dat$pork, dat$pork==(-10), 0)

# Cheese intake
dat$cheese <- dat$'V1408-0.0'
dat$cheese <- replace(dat$cheese, dat$cheese<0, NA)
dat$cheese <- replace(dat$cheese, dat$cheese==(-10), 0)

# Salt added to food
dat$salt <- dat$'V1478-0.0'
dat$salt <- replace(dat$salt, dat$salt<0, NA)
dat$salt <- replace(dat$salt, dat$salt==(-10), 0)

# Tea intake
dat$tea <- dat$'V1488-0.0'
dat$tea <- replace(dat$tea, dat$tea<0, NA)
dat$tea <- replace(dat$tea, dat$tea==(-10), 0)

# Coffee intake
dat$coffee <- dat$'V1498-0.0'
dat$coffee <- replace(dat$coffee, dat$coffee<0, NA)
dat$coffee <- replace(dat$coffee, dat$coffee==(-10), 0)

Efood <- c("veg_cook", "fish_oily", "fish_lean", "meat_proc", "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea", "coffee")
data.frame(Efood,Efood%in%colnames(dat))











##### Alcohol #####

# Alcohol intake frequency
dat$alc1 <- dat$'V1558-0.0'
dat$alc1 <- replace(dat$alc1, dat$alc1<0, NA)

Ealcohol <- c("alc1")
data.frame(Ealcohol,Ealcohol%in%colnames(dat))

# Final selection
allvars <- c(phen, Ivars, Ebody, Eill, Elife, Ephys, Esleep, Esmoke, Efood, Ealcohol)
dt <- dat[allvars]
dim(dt)

# Save output
saveRDS(dt, file = out_rdata)
cat("Saved:", out_rdata, "\n")
