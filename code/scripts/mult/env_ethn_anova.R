# .R----
.libPaths('/data2/morgante_lab/kgoda/software/R/x86_64-pc-linux-gnu-library/4.2')
.libPaths('/data2/morgante_lab/nklimko/software/R/x86_64-pc-linux-gnu-library/4.2')

library(data.table)
library(argparse)

set.seed(123)

mainCall <- function(dartpath,
                     envname,
                     outpath){
  # Demo data
  # dartpath <- 'data/mult/anova_input.rds'
  # envname <- 'Townsend'

  dart <- readRDS(dartpath)

  # Ensure ethnicity is a factor
  dart$Ethnicity <- as.factor(dart$Ethnicity)

  # Create formula for each environment
  freeform <- paste0(envname,' ~ Ethnicity')

  # Linear model
  tempo <- lm(data = dart, formula=freeform)

  # perform anova on model fit
  final <- anova(tempo)
  # final <- summary(tempo)

  # debug print
  print(tempo)

  # redundancy print
  print(final)

  # save to file
  saveRDS(final, file = outpath)

}

#args----
parser <- ArgumentParser(description= 'snakemake transfer')

parser$add_argument('--datapath')
parser$add_argument('--env')
parser$add_argument('--outpath')

snake <- parser$parse_args()
print(str(snake))

#call----
mainCall(snake$datapath,
         snake$env,
         snake$outpath)
