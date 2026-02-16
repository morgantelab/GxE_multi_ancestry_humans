### 2025/11/20
library(data.table)

# scrapes top enps for each environment based on specified cutoff (hline)
monet <- function(trait, env, hline){

  # ensure trait is allcaps for file system
  trait <- toupper(trait)

  # read in data using trait and environment combination
  newp <- readRDS(paste0('output/gpcxe/tables_gxegwas/', trait,'/',env,'.rds'))

  # select all snps above horizontal line determined from analysis/qqtop.Rmd qqplot analysis
  snpcheck <- newp[neglogP > hline,]

  # return data table of top snps
  final <- data.table(trait=trait, env=env, snpcheck)

  return(final)
}

## Diastolic
d1 <- monet('DP', 'smoking_now', 5.5)
d2 <- monet('DP', 'poultry', 5.2)
d3 <- monet('DP', 'cheese', 4.14)
d4 <- monet('DP', 'tea', 4.9)
d5 <- monet('DP', 'BFP', 4.7)

## Systolic
s1 <- monet('SP', 'TVtime', 5.5)
s2 <- monet('SP', 'sleep_d', 3.94)
s3 <- monet('SP', 'veg_cook', 4.7)
s4 <- monet('SP', 'pork', 4.23)
s5 <- monet('SP', 'tea', 4.155)
s6 <- monet('SP', 'coffee', 4.7)
s7 <- monet('SP', 'sleep_dev', 4.9)

## Pulse
p1 <- monet('PP', 'TVtime', 4.33)
p2 <- monet('PP', 'sleep_dev', 3.5)
p3 <- monet('PP', 'salt', 4.5)
p4 <- monet('PP', 'tea', 5.5)

# code to write out d1...p4 shown below
fim <- c(paste0('d',1:4), paste0('s',1:7), paste0('p',1:4))

fink <- paste0(fim, collapse=",")

# merge all outputs together
grands <- data.table(rbind(d1,d2,d3,d4,s1,s2,s3,s4,s5,s6,s7,p1,p2,p3,p4))

# write to file
fwrite(grands, 'output/gpcxe/topsnps.txt', sep='\t')

