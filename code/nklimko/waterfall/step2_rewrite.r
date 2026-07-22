set.seed(1123)
rm(list=ls())


library(optparse)
library(data.table)
library(dplyr)

# --- Parse options ---
option_list <- list(
  make_option("--input_rds", type = "character"),
  make_option("--withdrawn", type = "character"),
  make_option("--missing_geno", type = "character"),
  make_option("--output_rds", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$input_rds <- '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data1_20250501_only_analysis_vars.rds'
opt$withdrawn <- "/data2/morgante_lab/data/ukbiobank/ind_to_remove/withdraw62347_281_20240209.txt"
opt$missing_geno <- '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/missing_ids.txt'

# --- Load inputs ---
dtt <- readRDS(opt$input_rds)
withdrawn_ids <- fread(opt$withdrawn, header = FALSE)$V1
missing_ids <- fread(opt$missing_geno, header = FALSE)$V1


dtt <- dtt %>%
  mutate(
    ethn1_consolidated = case_when(
      ethn1_white == 1 ~ "White",
      ethn1_black == 1 ~ "Black",
      ethn1_chinese == 1 ~ "Chinese",
      ethn1_asian == 1 ~ "Asian",
      ethn1_mixed == 1 ~ "Mixed",
      TRUE ~ 'zz_other'  # if no value matches, set it as NA
    )
  )
r2 <- data.table(table(dtt$ethn1_consolidated))
i <- 2
names(r2)[2] <- 'base'

nrow(dtt)

# --- Subset by ID ---
dtt <- subset(dtt, !(ID %in% withdrawn_ids))
dtt <- subset(dtt, !(ID %in% missing_ids))


r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'genotyped'

nrow(dtt)



dtt$add.med10 <- dtt$add.med15 <- 0
dtt$add.med10 <- replace(dtt$add.med10, dtt$mediblood1==1 | dtt$mediblood2==1 | dtt$mediblood3==1 | dtt$mediblood4==1, 10)
dtt$add.med15 <- replace(dtt$add.med15, dtt$mediblood1==1 | dtt$mediblood2==1 | dtt$mediblood3==1 | dtt$mediblood4==1, 15)
dtt$SP0a <- dtt$SP0a+dtt$add.med15
dtt$DP0a <- dtt$DP0a+dtt$add.med10
dtt$PP0a <- dtt$SP0a-dtt$DP0a




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

dtt2 <- dtt
dtt <- dtt2


dtt <- subset(dtt, !is.na(DP0)); nrow(dtt) # diastolic blood pressure, first reading
dtt <- subset(dtt, !is.na(SP0)); nrow(dtt) # sistolic blood pressure, first reading
dtt <- subset(dtt, !is.na(PP0)); nrow(dtt) # pulse pressure, first reading

r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'DP_SP'

### Subset for phenotypes ###
#Removing data away by more than 4 sd for SP and DP
dtt <- dtt[abs(dtt$SP0 - mean(dtt$SP0, na.rm = TRUE)) <= 4 * sd(dtt$SP0, na.rm = TRUE), ]
dtt <- dtt[abs(dtt$DP0 - mean(dtt$DP0, na.rm = TRUE)) <= 4 * sd(dtt$DP0, na.rm = TRUE), ]

r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'DP_SP_4sd'


# Subset to remove rows with all mediblood columns equal to 0
dtt <- dtt[!(dtt$mediblood1 == 0 & dtt$mediblood2 == 0 & dtt$mediblood3 == 0 &
               dtt$mediblood4 == 0 & dtt$mediblood5 == 0 & dtt$mediblood6 == 0 &
               dtt$mediblood7 == 0), ]

r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'med_history'

### Subsetting for I vars ###
dtt <- subset(dtt, !is.na(Sex_SI)); nrow(dtt) # sex, self identified
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'self_sex'
dtt <- subset(dtt, !is.na(AOP)); nrow(dtt) # age at phenotyping
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'age_pheno'
dtt <- subset(dtt, !is.na(AOR)); nrow(dtt) # age at recruiment
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'age_recruit'
# Elife
dtt <- subset(dtt, !is.na(Townsend)); nrow(dtt) # Townsend deprivation index at recruitment
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'townsend'
# Ephys
dtt <- subset(dtt, !is.na(act0_d)); nrow(dtt) # number of days/week of moderate physical activity 10þ minutes
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'act0_d'
dtt <- subset(dtt, !is.na(TVtime)); nrow(dtt) # time spent watching television
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'TVtime'

# Esleep
dtt <- subset(dtt, !is.na(sleep_d)); nrow(dtt) # sleep duration
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'sleep'

# Esmoke
dtt <- subset(dtt, !is.na(smoking_now)); nrow(dtt) # is current smoker - smoking status
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'smoke_current'

# Efood
dtt <- subset(dtt, !is.na(veg_cook)); nrow(dtt) # cooked vegetable intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'veg_cook'

dtt <- subset(dtt, !is.na(fish_oily)); nrow(dtt) # oily fish intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'fish_oily'
dtt <- subset(dtt, !is.na(fish_lean)); nrow(dtt) # non-oily fish intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'fish_lean'
dtt <- subset(dtt, !is.na(meat_proc)); nrow(dtt) # processed meat intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'meat_proc'
dtt <- subset(dtt, !is.na(poultry)); nrow(dtt) # poultry intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'poultry'
dtt <- subset(dtt, !is.na(beef)); nrow(dtt) # beef intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'beef'
dtt <- subset(dtt, !is.na(lamb)); nrow(dtt) # lamb intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'lamb'
dtt <- subset(dtt, !is.na(pork)); nrow(dtt) # pork intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'pork'
dtt <- subset(dtt, !is.na(cheese)); nrow(dtt) # cheese intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'cheese'
dtt <- subset(dtt, !is.na(salt)); nrow(dtt) # salt added to food
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'salt'
dtt <- subset(dtt, !is.na(tea)); nrow(dtt) # tea intake
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'tea'
# Ealcohol
dtt <- subset(dtt, !is.na(alc1)); nrow(dtt) # alcohol intake frequency
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'alcohol'

# Ebody
dtt <- subset(dtt, !is.na(waist)); nrow(dtt) # waist circumference
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'waist'
### Adding pressure if medication ###

nrow(dtt)

# Subset to remove rows with any NA in mediblood columns
dtt <- dtt[!is.na(dtt$mediblood1) & !is.na(dtt$mediblood2) & !is.na(dtt$mediblood3) &
             !is.na(dtt$mediblood4) & !is.na(dtt$mediblood5) & !is.na(dtt$mediblood6) &
             !is.na(dtt$mediblood7), ]
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'mediblood_again'
nrow(dtt)

### Additional variables ###

dtt <- subset(dtt, !is.na(getup)); nrow(dtt) ## getup - ease of getting up in the morning ##
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'getup'

dtt <- subset(dtt, !is.na(coffee)); nrow(dtt) ## coffee - coffee consumption ##
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'coffee'
dtt <- subset(dtt, !is.na(smoked_past)); nrow(dtt) ## smoked_past - smoked in the past ##
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'smoke_past'
dtt <- subset(dtt, !is.na(BFP)); nrow(dtt) ## BFP - body fat percentage ##
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'bfp'


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
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'age_bound'
## Creating a cohort indicator, pasting Sex_SI and AOPc2 ##
dtt$coh <- paste(dtt$Sex_SI, dtt$AOPc2, sep='_')
table(dtt$coh)

table(dtt$ethn1, useNA = "always")
dtt <- subset(dtt, !is.na(ethn1)); nrow(dtt)
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'ethn_na_removed'

table(dtt$ethn1_white, useNA = "always")
table(dtt$ethn1_whbri, useNA = "always")
table(dtt$ethn1_black, useNA = "always")
table(dtt$ethn1_asian, useNA = "always")
table(dtt$ethn1_mixed, useNA = "always")
table(dtt$ethn1_chinese, useNA = "always")
table(dtt$ethn1_other, useNA = "always")

# Subset to remove rows with other ethnicities.
dtt <- dtt[!(dtt$ethn1 == 6), ]
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'ethn_other_removed'
nrow(dtt)

#Removing 99% percentile for coffee, tea, veg_cook, and TVtime
dtt <- dtt[dtt$coffee <= quantile(dtt$coffee, 0.99), ]
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'coffee_99'
dtt <- dtt[dtt$tea <= quantile(dtt$tea, 0.99), ]
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'tea_99'
dtt <- dtt[dtt$veg_cook <= quantile(dtt$veg_cook, 0.99), ]
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'veg_cook_99'
dtt <- dtt[dtt$TVtime <= quantile(dtt$TVtime, 0.99), ]
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'tvtime_99'
#Removing 1% percentile and 99% percentile for sleep duration
dtt <- dtt[dtt$sleep_d >= quantile(dtt$sleep_d, 0.01), ]
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'sleep_01'
dtt <- dtt[dtt$sleep_d <= quantile(dtt$sleep_d, 0.99), ]
r2 <- data.table(table(dtt$ethn1_consolidated))[r2, on=.(V1), all=T]
i <- i + 1
names(r2)[2] <- 'sleep_99'

dtt <- dtt %>%
  mutate(
    ethn1_consolidated = case_when(
      ethn1_white == 1 ~ "White",
      ethn1_black == 1 ~ "Black",
      ethn1_chinese == 1 ~ "Chinese",
      ethn1_asian == 1 ~ "Asian",
      ethn1_mixed == 1 ~ "Mixed",
      TRUE ~ NA_character_  # if no value matches, set it as NA
    )
  )

dtt$sleep_dev <- (dtt$sleep - mean(dtt$sleep))^2



#
#
#    trait    env   CHR        BP         ID      BETA       SE           P  neglogP
#    <char> <char> <int>     <int>     <char>     <num>    <num>      <char>    <num>
# 1:     DP    BFP     1  12156236 rs11121855  0.651747 0.144619  6.6166e-06 5.179365
# 2:     DP    BFP     3 174348263   rs556852  0.494468 0.112534 1.11809e-05 4.951523
# 3:     DP    BFP     4   7985830  rs7686574 -0.716256 0.164630 1.36256e-05 4.865644
# 4:     DP    BFP     5 162213246 rs62398131  0.500224 0.110624 6.15969e-06 5.210441
# 5:     DP    BFP     6  89444829  rs9344849 -0.534620 0.119392 7.57459e-06 5.120641
# 6:     DP    BFP     6  89561806 rs34917906 -0.604514 0.132833 5.36647e-06 5.270311
# 7:     DP    BFP    10 131684974  rs7911984 -0.588792 0.123615 1.91798e-06 5.717156
# 8:     DP    BFP    14  20933806 rs77502768  0.712190 0.158567 7.10759e-06 5.148278
#
# g2 <- grands[env=='BFP',]
# g2$CHR <- NULL
# g2$BP <- NULL
# g2$BETA <- NULL
# g2$SE <- NULL
# g2$neglogP <- NULL
#
# g2$P <- as.numeric(g2$P)
# g2 <- g2[order(P),]
# g2$P <- scientific(g2$P, digits=3)
# g2
#
# DP & BFP & rs7911984 & 1.92e-06 \\
# DP & BFP & rs34917906 & 5.37e-06 \\
# DP & BFP & rs62398131 & 6.16e-06 \\
# DP & BFP & rs11121855 & 6.62e-06 \\
# DP & BFP & rs77502768 & 7.11e-06 \\
# DP & BFP & rs9344849 & 7.57e-06 \\
# DP & BFP & rs556852 & 1.12e-05 \\
# DP & BFP & rs7686574 & 1.36e-05 \\
