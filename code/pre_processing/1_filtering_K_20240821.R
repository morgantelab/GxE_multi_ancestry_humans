rm(list=ls()); gc()

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data")

library(data.table)

##### Load individuals to remove #####
inds_to_remove <- fread("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry_backup/march_run_updated_dataset/1_pre_processing_dataset/inputs/given_dataset_n_withdrew/given_dataset_FT/inds_to_withdraw_62347_281_20240209_FM.txt", data.table = FALSE)$V1
#length(inds_to_remove)
#head(inds_to_remove)

##### geno missin individuals to remove #######
inds_to_remove_geno <- fread("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry_backup/march_run_updated_dataset/1_pre_processing_dataset/inputs/missing_geno/missing_ids.txt", data.table = FALSE)$V1
#length(inds_to_remove_geno)
#head(inds_to_remove_geno)

##### Load variable assignement #####
load("data1_20240821.RData")
dtt <- dt
nrow(dtt)

### Subset by ID ###
dtt <- subset(dtt, !(ID%in%inds_to_remove))
nrow(dtt)
dtt <- subset(dtt, !(ID%in%inds_to_remove_geno))
nrow(dtt)

sum(is.na(dtt$PR0))
sum(is.na(dtt$DP0))
sum(is.na(dtt$SP0))
sum(is.na(dtt$PP0))
sum(is.na(dtt$DP0a))
sum(is.na(dtt$SP0a))

sum(is.na(dtt$Sex_SI))
sum(is.na(dtt$AOP))
sum(is.na(dtt$AOR))

sum(is.na(dtt$Townsend))
sum(is.na(dtt$act0_d))
sum(is.na(dtt$TVtime))
sum(is.na(dtt$sleep_d))
sum(is.na(dtt$smoking_now))
sum(is.na(dtt$veg_cook))
sum(is.na(dtt$fish_oily))
sum(is.na(dtt$fish_lean))
sum(is.na(dtt$meat_proc))
sum(is.na(dtt$poultry))
sum(is.na(dtt$beef))
sum(is.na(dtt$lamb))
sum(is.na(dtt$pork))
sum(is.na(dtt$cheese))
sum(is.na(dtt$salt))
sum(is.na(dtt$tea))
sum(is.na(dtt$alc1))
sum(is.na(dtt$waist))
sum(is.na(dtt$getup))
sum(is.na(dtt$coffee))
sum(is.na(dtt$smoked_past))
sum(is.na(dtt$BFP))


### Subset for phenotypes ###
#Removing data away by more than 4 sd for SP and DP
dtt <- dtt[abs(dtt$SP0 - mean(dtt$SP0, na.rm = TRUE)) <= 4 * sd(dtt$SP0, na.rm = TRUE), ]
dtt <- dtt[abs(dtt$DP0 - mean(dtt$DP0, na.rm = TRUE)) <= 4 * sd(dtt$DP0, na.rm = TRUE), ]

dtt <- subset(dtt, !is.na(PR0)); nrow(dtt) # pulse rate, automatic
dtt <- subset(dtt, !is.na(DP0)); nrow(dtt) # diastolic blood pressure, first reading
dtt <- subset(dtt, !is.na(SP0)); nrow(dtt) # sistolic blood pressure, first reading
dtt <- subset(dtt, !is.na(PP0)); nrow(dtt) # pulse pressure, first reading

### Subsetting for I vars ###
dtt <- subset(dtt, !is.na(Sex_SI)); nrow(dtt) # sex, self identified
dtt <- subset(dtt, !is.na(AOP)); nrow(dtt) # age at phenotyping
dtt <- subset(dtt, !is.na(AOR)); nrow(dtt) # age at recruiment

# Elife
dtt <- subset(dtt, !is.na(Townsend)); nrow(dtt) # Townsend deprivation index at recruitment

# Ephys
dtt <- subset(dtt, !is.na(act0_d)); nrow(dtt) # number of days/week of moderate physical activity 10þ minutes
dtt <- subset(dtt, !is.na(TVtime)); nrow(dtt) # time spent watching television

# Esleep
dtt <- subset(dtt, !is.na(sleep_d)); nrow(dtt) # sleep duration

# Esmoke
dtt <- subset(dtt, !is.na(smoking_now)); nrow(dtt) # is current smoker - smoking status

# Efood
dtt <- subset(dtt, !is.na(veg_cook)); nrow(dtt) # cooked vegetable intake
dtt <- subset(dtt, !is.na(fish_oily)); nrow(dtt) # oily fish intake
dtt <- subset(dtt, !is.na(fish_lean)); nrow(dtt) # non-oily fish intake
dtt <- subset(dtt, !is.na(meat_proc)); nrow(dtt) # processed meat intake
dtt <- subset(dtt, !is.na(poultry)); nrow(dtt) # poultry intake
dtt <- subset(dtt, !is.na(beef)); nrow(dtt) # beef intake
dtt <- subset(dtt, !is.na(lamb)); nrow(dtt) # lamb intake
dtt <- subset(dtt, !is.na(pork)); nrow(dtt) # pork intake
dtt <- subset(dtt, !is.na(cheese)); nrow(dtt) # cheese intake
dtt <- subset(dtt, !is.na(salt)); nrow(dtt) # salt added to food
dtt <- subset(dtt, !is.na(tea)); nrow(dtt) # tea intake

# Ealcohol
dtt <- subset(dtt, !is.na(alc1)); nrow(dtt) # alcohol intake frequency

# Ebody
dtt <- subset(dtt, !is.na(waist)); nrow(dtt) # waist circumference

### Adding pressure if medication ###

dtt$add.med10 <- dtt$add.med15 <- 0
dtt$add.med10 <- replace(dtt$add.med10, dtt$mediblood1==1 | dtt$mediblood2==1 | dtt$mediblood3==1 | dtt$mediblood4==1, 10)
dtt$add.med15 <- replace(dtt$add.med15, dtt$mediblood1==1 | dtt$mediblood2==1 | dtt$mediblood3==1 | dtt$mediblood4==1, 15)
dtt$SP0a <- dtt$SP0a+dtt$add.med15
dtt$DP0a <- dtt$DP0a+dtt$add.med10
dtt$PP0a <- dtt$SP0a-dtt$DP0a

# Subset to remove rows with all mediblood columns equal to 0
dtt <- dtt[!(dtt$mediblood1 == 0 & dtt$mediblood2 == 0 & dtt$mediblood3 == 0 & 
               dtt$mediblood4 == 0 & dtt$mediblood5 == 0 & dtt$mediblood6 == 0 & 
               dtt$mediblood7 == 0), ]
nrow(dtt)

# Subset to remove rows with any NA in mediblood columns
dtt <- dtt[!is.na(dtt$mediblood1) & !is.na(dtt$mediblood2) & !is.na(dtt$mediblood3) & 
             !is.na(dtt$mediblood4) & !is.na(dtt$mediblood5) & !is.na(dtt$mediblood6) & 
             !is.na(dtt$mediblood7), ]

nrow(dtt)

### Additional variables ###

dtt <- subset(dtt, !is.na(getup)); nrow(dtt) ## getup - ease of getting up in the morning ##
dtt <- subset(dtt, !is.na(coffee)); nrow(dtt) ## coffee - coffee consumption ##
dtt <- subset(dtt, !is.na(smoked_past)); nrow(dtt) ## smoked_past - smoked in the past ##
dtt <- subset(dtt, !is.na(BFP)); nrow(dtt) ## BFP - body fat percentage ##



###################################
###  Creating design indicators ###
###################################

### Age at phenotyping classes: AOPc2 ###

## Boundary values taken and adjusted from Kerin and Marchini, 2020 ##
summary(dtt$AOP)
x1 <- 40; x1
x2 <- 51; x2
x3 <- 58; x3
x4 <- 63; x4
x5 <- 70; x5
dtt$AOPc2 <- 0
dtt$AOPc2[which(dtt$AOP>=x1)] <- 1
dtt$AOPc2[which(dtt$AOP>=x2)] <- 2
dtt$AOPc2[which(dtt$AOP>=x3)] <- 3
dtt$AOPc2[which(dtt$AOP>=x4)] <- 4
dtt$AOPc2[which(dtt$AOP>x5)] <- 5
table(dtt$AOPc2)

## Keeping only between 40 and 70 years old at phenotyping ##
dtt <- subset(dtt, AOPc2!=0); nrow(dtt)
dtt <- subset(dtt, AOPc2!=5); nrow(dtt)

## Creating a cohort indicator, pasting Sex_SI and AOPc2 ##
dtt$coh <- paste(dtt$Sex_SI, dtt$AOPc2, sep='_')
table(dtt$coh)

table(dtt$ethn1, useNA = "always")
dtt <- subset(dtt, !is.na(ethn1)); nrow(dtt)

table(dtt$ethn1_white, useNA = "always")
table(dtt$ethn1_whbri, useNA = "always")
table(dtt$ethn1_black, useNA = "always")
table(dtt$ethn1_asian, useNA = "always")
table(dtt$ethn1_mixed, useNA = "always")
table(dtt$ethn1_chinese, useNA = "always")
table(dtt$ethn1_other, useNA = "always")

# Subset to remove rows with other ethnicities.
dtt <- dtt[!(dtt$ethn1 == 6), ]
nrow(dtt)

#Removing 99% percentile for coffee, tea, veg_cook, and TVtime
dtt <- dtt[dtt$coffee <= quantile(dtt$coffee, 0.99), ]
dtt <- dtt[dtt$tea <= quantile(dtt$tea, 0.99), ]
dtt <- dtt[dtt$veg_cook <= quantile(dtt$veg_cook, 0.99), ]
dtt <- dtt[dtt$TVtime <= quantile(dtt$TVtime, 0.99), ]

#Removing 1% percentile and 99% percentile for sleep duration
dtt <- dtt[dtt$sleep_d >= quantile(dtt$sleep_d, 0.01), ]
dtt <- dtt[dtt$sleep_d <= quantile(dtt$sleep_d, 0.99), ]

#Creating a new variable called sleep_sd
#to measure how far is sleep duration value from population mean sleep duration
#dtt$sleep_sd <- abs((dtt$sleep_d - mean(dtt$sleep_d, na.rm = TRUE)) / sd(dtt$sleep_d, na.rm = TRUE))
#sleep_sd <- data.table(scale((data$sleep - mean(data$sleep))^2))
#dtt$sleep_sd <- 

## Creating 8 random subsets, same number as the cohorts ##
dtt$sub <- sample(11:18, nrow(dtt), replace=T)
table(dtt$sub)

## Checking the distribution of the random samples ##
table(dtt$sub, dtt$Sex_SI)
table(dtt$sub, dtt$AOPc2)
table(dtt$sub, dtt$coh)

#below doesnt work
 ## Checking the distribution of the Tiezzi-approved samples ##
 # > table(dtt$tiezzi, dtt$Sex_SI)
#Error in table(dtt$tiezzi, dtt$Sex_SI) : 
#  all arguments must have the same length
#> table(dtt$tiezzi, dtt$AOPc2)
#Error in table(dtt$tiezzi, dtt$AOPc2) : 
#  all arguments must have the same length
#> table(dtt$tiezzi, dtt$coh)
#Error in table(dtt$tiezzi, dtt$coh) : 
#  all arguments must have the same length

## Checking the distribution of the Tiezzi-approved samples ##
#table(dtt$tiezzi, dtt$Sex_SI)
#table(dtt$tiezzi, dtt$AOPc2)
#table(dtt$tiezzi, dtt$coh)

save(dtt, file="data2_20240930.RData")

rm(list=ls()); gc()