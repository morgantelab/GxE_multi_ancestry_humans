
library(data.table)
library(TeachingDemos)
# rebirth ----

trait <- 'DP'
model <- 'Model 6a'
datapath <-  'output/model6_alt/DP/pcrelate/varcomp.csv'

trait <- 'DP'
model <- 'Model 1'
datapath <- 'data/model/varabs_DP_pcrelate_just_X1.csv'

wreath <- function(veclist){
  trait <- veclist[1]
  model <- veclist[2]
  datapath <- veclist[3]

  temp2 <- fread(datapath)
  temp2$V1 <- NULL

  names(temp2) <- gsub('V_','',names(temp2))
  temp3 <- suppressWarnings(melt(temp2))

  stats8 <- temp3 %>%
    group_by(variable) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      lower = emp.hpd(value, conf=0.95)[1],
      upper = emp.hpd(value, conf=0.95)[2],
      .groups = "drop"
    )

  final <- data.table(stats8)

  names(final)[1] <- 'comp'
  final$trait <- trait
  final$model <- model
  final$modnum <- substr(model,7,7)

  return(final)
}

palorse <- function(veclist){
  trait <- veclist[1]
  modnum <- veclist[2]
  datapath <- veclist[3]

  temp2 <- fread(datapath)
  temp2$V1 <- NULL

  names(temp2) <- gsub('V_', paste0(trait,'_',modnum,'_'), names(temp2))

  return(temp2)

}

temp2

mcount <- 8
traitlist <-
  c(rep('DP', mcount),
    rep('SP', mcount),
    rep('PP', mcount))

numlist <- rep(c(1:7, '6b'), 3)

filelist <- c('data/model/varabs_DP_pcrelate_just_X1.csv',
              'data/model/varabs_DP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/DP/pcrelate/varcomp.csv',
              'output/brit_24k/DP/plink/model6_white24k.csv',
              'data/model/varabs_SP_pcrelate_just_X1.csv',
              'data/model/varabs_SP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/SP/pcrelate/varcomp.csv',
              'output/brit_24k/SP/plink/model6_white24k.csv',
              'data/model/varabs_PP_pcrelate_just_X1.csv',
              'data/model/varabs_PP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/PP/pcrelate/varcomp.csv',
              'output/brit_24k/PP/plink/model6_white24k.csv')


hopelist <- data.table(traitlist, modlist, filelist)
hopelist <- data.table(traitlist, numlist, filelist)


gg <- vector(mode='list', length=nrow(hopelist))

for( i in 1:nrow(hopelist)){
  gg[[i]] <- palorse(unlist(hopelist[i,1:3,with=0]))
}

# cat(paste0('gg[[',1:24,']], '))
velet <- cbind(gg[[1]],  gg[[2]],  gg[[3]],  gg[[4]],  gg[[5]],  gg[[6]],  gg[[7]],  gg[[8]],  gg[[9]],  gg[[10]],  gg[[11]],  gg[[12]],  gg[[13]],  gg[[14]],  gg[[15]],  gg[[16]],  gg[[17]],  gg[[18]],  gg[[19]],  gg[[20]],  gg[[21]],  gg[[22]],  gg[[23]],  gg[[24]])

malot <- ess(velet)


melt(malot)
midas <- data.table(names(velet), malot)
names(midas) <- c('chord', 'ess')


midas$ess_int <- round(midas$ess, 0)
midas$ess_int <- NULL
midas$trait <- str_split_i(midas$chord, '_', 1)
midas$model <- str_split_i(midas$chord, '_', 2)
midas$component <- str_split_i(midas$chord, '_', 3)
midas$chord <- NULL
setcolorder(midas, c(2:4, 1))
kable(midas, digits=0, 'simple')

demo <- rbindlist(gg)
esso <- cbindlist(gg)
fwrite(demo, 'data/nklimko/demo_varcomp_1-7.txt', sep='\t')

mean(midas$ess)
midas[component=='GxE', mean(ess)]
midas[component=='GxE', mean(ess), by=.(trait)]

# demo <- wreath('DP', 'Model 6a', 'output/model6_alt/DP/pcrelate/varcomp.csv')

gg <- vector(mode='list', length=nrow(hopelist))

for( i in 1:nrow(hopelist)){
  gg[[i]] <- wreath(unlist(hopelist[i,1:3,with=0]))
}


demo <- rbindlist(gg)
fwrite(demo, 'data/nklimko/demo_varcomp_1-7.txt', sep='\t')

library(stringr)

# str_split_i
demo$modnum <- str_split_i(demo$model, ' ', 2)

demo$model <- gsub('Model 8', 'Model 6 White British',demo$model)
# demo$model2 <- NULL
demo <- demo[order(trait, model),]
sixes <- demo[modnum ==6 | modnum==8,]
sixes$modnum <- NULL

kable(sixes, 'simple')

fwrite(demo, 'data/dump/model_components.txt', sep='\t')

library(data.table)
dart <- fread('data/dump/model_components.txt')

dart[trait=='DP',]
dart[trait=='SP',]
dart[trait=='PP',]


# ashes ----

if(0){
  # setwd('')

  # existing <- read.csv("data/model/figure_4_dataset.csv")
  # existingold <- fread('data/model/vc_means_all_models_with_ci.csv')
  #
  # existing <- fread('data/model/vc_means_all_models_with_ci.csv')
  #
  #
  # datapath <- 'output/model7/PP/plink/varcomp.csv'
  # temp <- fread('output/model8/PP/plink/varcomp.csv')
  # temp2 <- fread('output/model8/PP/plink/varcomp.csv')

  library(data.table)
  wreath <- function(trait, model, datapath){
    temp2 <- fread(datapath)
    temp2$V1 <- NULL

    names(temp2) <- gsub('V_','',names(temp2))
    temp3 <- suppressWarnings(melt(temp2))

    stats8 <- temp3 %>%
      group_by(variable) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        lower = mean(value, na.rm = TRUE) - qt(0.975, sum(!is.na(value)) - 1) * sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
        upper = mean(value, na.rm = TRUE) + qt(0.975, sum(!is.na(value)) - 1) * sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
        .groups = "drop"
      )

    final <- data.table(stats8)

    names(final)[1] <- 'comp'
    final$trait <- trait
    final$model <- model

    return(final)
  }

}


# demon$mean <- round(demon$mean, 2)
# demon$lower <- round(demon$lower, 2)
# demon$upper <- round(demon$upper, 2)
#
# demon[order(trait, comp, model),]


mcount <- 8
mcount <- 7

traitlist <-
  c(rep('DP', mcount),
    rep('SP', mcount),
    rep('PP', mcount))


modlist <- rep(paste0('Model ',1:mcount), 3)
filelist <- c('data/model/varabs_DP_pcrelate_just_X1.csv',
              'data/model/varabs_DP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/DP/pcrelate/varcomp.csv',
              'data/model/varabs_SP_pcrelate_just_X1.csv',
              'data/model/varabs_SP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/SP/pcrelate/varcomp.csv',
              'data/model/varabs_PP_pcrelate_just_X1.csv',
              'data/model/varabs_PP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/PP/pcrelate/varcomp.csv')

filelist <- c('data/model/varabs_DP_pcrelate_just_X1.csv',
              'data/model/varabs_DP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/DP/pcrelate/varcomp.csv',
              'output/brit_24k/DP/plink/model6_white24k.csv',
              'data/model/varabs_SP_pcrelate_just_X1.csv',
              'data/model/varabs_SP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/SP/pcrelate/varcomp.csv',
              'output/brit_24k/SP/plink/model6_white24k.csv',
              'data/model/varabs_PP_pcrelate_just_X1.csv',
              'data/model/varabs_PP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/PP/pcrelate/varcomp.csv',
              'output/brit_24k/PP/plink/model6_white24k.csv')

filelist <- c('data/model/varabs_DP_pcrelate_just_X1.csv',
              'data/model/varabs_DP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_DP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/DP/pcrelate/varcomp.csv',
              'output/model8/DP/pcrelate/varcomp.csv',
              'output/model9/DP/pcrelate/varcomp.csv',
              'data/model/varabs_SP_pcrelate_just_X1.csv',
              'data/model/varabs_SP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_SP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model8/SP/pcrelate/varcomp.csv',
              'output/model7/SP/pcrelate/varcomp.csv',
              'output/model9/SP/pcrelate/varcomp.csv',
              'data/model/varabs_PP_pcrelate_just_X1.csv',
              'data/model/varabs_PP_pcrelate_just_X1_X2.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_E.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G_E.csv',
              'data/model/varabs_PP_pcrelate_X1_X2_G_E_GxE.csv',
              'output/model7/PP/pcrelate/varcomp.csv',
              'output/model8/PP/pcrelate/varcomp.csv',
              'output/model9/PP/pcrelate/varcomp.csv')


hopelist <- data.table(traitlist, modlist, filelist)
#
# d7 <- wreath('DP', 'Model 7', 'output/model7/DP/pcrelate/varcomp.csv')
# d8 <- wreath('DP', 'Model 8', 'output/model8/DP/pcrelate/varcomp.csv')
# d9 <- wreath('DP', 'Model 9', 'output/model9/DP/pcrelate/varcomp.csv')
#
# s8 <- wreath('SP', 'Model 8', 'output/model8/SP/pcrelate/varcomp.csv'
#              s7 <- wreath('SP', 'Model 7', 'output/model7/SP/pcrelate/varcomp.csv'
#                           s9 <- wreath('SP', 'Model 9', 'output/model9/SP/pcrelate/varcomp.csv'
#
#                                        p7 <- wreath('PP', 'Model 7', 'output/model7/PP/pcrelate/varcomp.csv')
#                                        p8 <- wreath('PP', 'Model 8', 'output/model8/PP/pcrelate/varcomp.csv')
#                                        p9 <- wreath('PP', 'Model 9', 'output/model9/PP/pcrelate/varcomp.csv')
#

if(1){



  data/model/varabs_PP_pcrelate_X1_X2_G_E_GxE.csv


  d6a <- wreath('DP', 'Model 6a', 'output/model6_alt/DP/pcrelate/varcomp.csv')
  s6a <- wreath('SP', 'Model 6a', 'output/model6_alt/SP/pcrelate/varcomp.csv')
  p6a <- wreath('PP', 'Model 6a', 'output/model6_alt/PP/pcrelate/varcomp.csv')
  d7 <- wreath('DP', 'Model 7', 'output/model7/DP/pcrelate/varcomp.csv')
  s7 <- wreath('SP', 'Model 7', 'output/model7/SP/pcrelate/varcomp.csv')
  p7 <- wreath('PP', 'Model 7', 'output/model7/PP/pcrelate/varcomp.csv')
  d8 <- wreath('DP', 'Model 8', 'output/model8/DP/pcrelate/varcomp.csv')
  s8 <- wreath('SP', 'Model 8', 'output/model8/SP/pcrelate/varcomp.csv')
  p8 <- wreath('PP', 'Model 8', 'output/model8/PP/pcrelate/varcomp.csv')
  d9 <- wreath('DP', 'Model 9', 'output/model9/DP/pcrelate/varcomp.csv')
  s9 <- wreath('SP', 'Model 9', 'output/model9/SP/pcrelate/varcomp.csv')
  p9 <- wreath('PP', 'Model 9', 'output/model9/PP/pcrelate/varcomp.csv')


  d6w <- wreath('DP', 'Model 6w', 'output/white_24k/DP/plink/model6_white24k.csv')
  s6w <- wreath('SP', 'Model 6w', 'output/white_24k/SP/plink/model6_white24k.csv')
  p6w <- wreath('PP', 'Model 6w', 'output/white_24k/PP/plink/model6_white24k.csv')

  d6b <- wreath('DP', 'Model 6b', 'output/brit_24k/DP/plink/model6_white24k.csv')
  s6b <- wreath('SP', 'Model 6b', 'output/brit_24k/SP/plink/model6_white24k.csv')
  p6b <- wreath('PP', 'Model 6b', 'output/brit_24k/PP/plink/model6_white24k.csv')

}


# d7 <- wreath('DP', 'Model 7', 'output/model7/DP/plink/varcomp.csv')
# s7 <- wreath('SP', 'Model 7', 'output/model7/SP/plink/varcomp.csv')
# p7 <- wreath('PP', 'Model 7', 'output/model7/PP/plink/varcomp.csv')
# d8 <- wreath('DP', 'Model 8', 'output/model8/DP/plink/varcomp.csv')
# s8 <- wreath('SP', 'Model 8', 'output/model8/SP/plink/varcomp.csv')
# p8 <- wreath('PP', 'Model 8', 'output/model8/PP/plink/varcomp.csv')

# if(0){
#   demon <- rbind(d6a, s6a, p6a, d7, s7, p7, d8, s8, p8)
# }
# demo <- wreath('PP', 'Model 8', 'output/model8/PP/plink/varcomp.csv')
# demo


# existing figure 4 dataset ----
if(1){
  existing <- read.csv("data/model/figure_4_dataset.csv")
  # existing$colname <- gsub('V_','',existing$colname)
  # existing$trait <- substr(existing$filename, 1,2)
  # existing$model <- substr(existing$filename, 4, nchar(existing$filename))
  # existing$filename <- NULL
  # names(existing) <- c('comp', 'mean', 'lower', 'upper', 'trait', 'model')

  names(existing) <- c('trait', 'model', 'comp', 'mean', 'lower', 'upper')

  existing$comp <- gsub('V_','',existing$comp)
  existing$model <- gsub('V_','',existing$model)

  setcolorder(existing, c(3:6,1:2))
  unique(existing$model)

  existing$model <- gsub("^X1$", 'Model 1', existing$model)
  existing$model <- gsub("^X1_X2$", 'Model 2', existing$model)
  existing$model <- gsub("^X1_X2_G$", 'Model 3', existing$model)
  existing$model <- gsub("^X1_X2_E$", 'Model 4', existing$model)
  existing$model <- gsub("^X1_X2_G_E$", 'Model 5', existing$model)
  existing$model <- gsub("^X1_X2_G_E_GxE$", 'Model 6', existing$model)
  # unique(existing$model)
}

# subject matter----
if(1){
  matte <- demo

  matte <- rbind(existing, d6a, s6a, p6a, d6w, s6w, p6w, d6b, s6b, p6b, d7, s7, p7, d8, s8, p8, d9, s9, p9)

  matte$model <- gsub('GxE', 'Model 6', matte$model)
  matte$model <- gsub('G_E', 'Model 5', matte$model)
  matte$model <- gsub('E', 'Model 4', matte$model)
  matte$model <- gsub('G', 'Model 3', matte$model)
  matte$model <- gsub('X1_X2', 'Model 2', matte$model)
  matte$model <- gsub('X1', 'Model 1', matte$model)

  matte$comp2 <- matte$comp

  matte$comp2 <- gsub('^X1$','A_X1',matte$comp2)
  matte$comp2 <- gsub('^X2$','B_X2',matte$comp2)
  matte$comp2 <- gsub('^E$','C_E',matte$comp2)
  matte$comp2 <- gsub('^G$','D_G',matte$comp2)
  matte$comp2 <- gsub('^m7$','E_M7',matte$comp2)
  matte$comp2 <- gsub('^GxE$','F_GxE',matte$comp2)

  # matte$fax <- paste0(matte$comp,'_from_',matte$model)
  matte$fax <- paste0(matte$model, '_', matte$comp2)
  matte$modnum <- gsub('Model ', '', matte$model)

  matte <- data.table(matte)

  mat_dp <- matte[trait=='DP',]
  mat_sp <- matte[trait=='SP',]
  mat_pp <- matte[trait=='PP',]
}

# fwrite(matte, file='data/nklimko/matte.txt')

# table(matte$comp2)

# values <- c('A_X1', 'B_X2', 'C_E', 'D_G', 'E_M7', 'F_GxE'), labels <- c('X1', 'X2', 'E', 'G', 'M7', 'GxE')
# subject matter function----
matter <- function(datanow, fulltitle, montag=NULL, mod7=TRUE){
  plothole <- ggplot(data=datanow, aes(x=fax, y=mean))+
    # scale_color_manual(discrete=T, values = c('A_X1', 'B_X2', 'C_E', 'D_G', 'E_M7', 'F_GxE'), labels = c('X1', 'X2', 'E', 'G', 'M7', 'GxE'))+
    {if(mod7)  scale_color_discrete(labels = c('X1', 'X2', 'E', 'G', 'M7', 'GxE'))}+
    {if(!mod7) scale_color_discrete(labels = c('X1', 'X2', 'E', 'G', 'GxE'))}+
    geom_point(aes(color=comp2))+
    # geom_text(aes(label=modnum))+
    geom_errorbar(aes(ymin = lower, ymax = upper, color=comp2),
                  width = 0.8, linewidth = 0.4, alpha = 1) +
    # labs(title='Model 8 Variance Components')+
    labs(title=fulltitle, color='Component', tag=montag)+
    xlab('Model')+
    # xlab(NULL)+
    ylab('PoVE')+
    theme_minimal()+
    theme(axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())+
    geom_vline(xintercept = verty, linetype=2, color='lightgray')+
    geom_text(data=modcode, aes(label=mn, x = x, y = y), inherit.aes=FALSE)

  return(plothole)
}

#plot setup----
if(1){
  # all of this sets up vertical dashed lines to split models manually

  #number of models
  # modelnumber <- 1:length(unique(matte$model))

  # indexes of line breaks
  # verty <- c(1.5, 3.5, 6.5, 9.5, 13.5, 17.5, 22.5, 27.5, 33.5)
  verty <- 0.5 + cumsum(unlist(table(mat_dp$model)))
  verty <- verty[-length(verty)]

  # add on hardpoint ends
  verty2 <- c(0, verty, nrow(mat_dp)+1)

  # calculates midpoint of each section to place Model Number below x axis
  verty3 <- (verty2[2:length(verty2)] - verty2[1:(length(verty2) - 1)]) / 2 + verty2[1:(length(verty2) - 1)]


  # altnumber <- c(1:6, '6a', '6b', '6w', 7:9)
  altnumber <- 1:mcount

  #matter function plots model number at specified breakpoint at y = -1 for each label.
  modcode <- data.table(mn=altnumber, x=verty3, y=-1)

}


#results maker----
if(1){
  mat_dp <- matte[trait=='DP',]
  mat_sp <- matte[trait=='SP',]
  mat_pp <- matte[trait=='PP',]

  dop <- matter(mat_dp, fulltitle=NULL, montag='DP')
  sop <- matter(mat_sp, fulltitle=NULL, montag='SP')
  pop <- matter(mat_pp, fulltitle=NULL, montag='PP')
  # dop <- matter(mat_dp, fulltitle='Diastolic Pressure')
  # sop <- matter(mat_sp, fulltitle='Systolic Pressure')
  # pop <- matter(mat_pp, fulltitle='Pulse Pressure')
  # print(plot_grid(dop, sop, pop, ncol=1))
}

print(plot_grid(dop, sop, pop, ncol=1))

# labels <- c('A_X1' = 'X1',
# 'B_X2' = 'X2',
# 'C_E' = 'E',
# 'D_G' = 'G',
# 'E_M7' = 'M7',
# 'F_GxE' = 'GxE')



# fin_dp


# modcode

if(0){

  # who up profounding they hearing loss


  mattix <- mat_pp[modnum==6 | modnum ==8,]




  # Clean environment
  rm(list = ls())
  set.seed(1123)

  .libPaths('/data2/morgante_lab/kgoda/software/R/x86_64-pc-linux-gnu-library/4.2')

  library(dplyr)
  library(tidyr)

  setwd("")

  # Define the list of files to process
  file_paths <- list(
    DP_X1 = "data/model/varabs_DP_pcrelate_just_X1.csv",
    SP_X1 = "data/model/varabs_SP_pcrelate_just_X1.csv",
    PP_X1 = "data/model/varabs_PP_pcrelate_just_X1.csv",
    DP_X1_X2 = "data/model/varabs_DP_pcrelate_just_X1_X2.csv",
    SP_X1_X2 = "data/model/varabs_SP_pcrelate_just_X1_X2.csv",
    PP_X1_X2 = "data/model/varabs_PP_pcrelate_just_X1_X2.csv",
    DP_G_E = "data/model/varabs_DP_pcrelate_X1_X2_G_E.csv",
    SP_G_E = "data/model/varabs_SP_pcrelate_X1_X2_G_E.csv",
    PP_G_E = "data/model/varabs_PP_pcrelate_X1_X2_G_E.csv",
    DP_G = "data/model/varabs_DP_pcrelate_X1_X2_G.csv",
    SP_G = "data/model/varabs_SP_pcrelate_X1_X2_G.csv",
    PP_G = "data/model/varabs_PP_pcrelate_X1_X2_G.csv",
    DP_GxE = "data/model/varabs_DP_pcrelate_X1_X2_G_E_GE.csv",
    SP_GxE = "data/model/varabs_SP_pcrelate_X1_X2_G_E_GE.csv",
    PP_GxE = "data/model/varabs_PP_pcrelate_X1_X2_G_E_GE.csv",
    DP_E = "data/model/varabs_DP_pcrelate_X1_X2_E.csv",
    SP_E = "data/model/varabs_SP_pcrelate_X1_X2_E.csv",
    PP_E = "data/model/varabs_PP_pcrelate_X1_X2_E.csv",
    DP7 = 'output/model7/DP/plink/varcomp.csv',
    SP7 = 'output/model7/SP/plink/varcomp.csv',
    PP7 = 'output/model7/PP/plink/varcomp.csv',
    DP8 = 'output/model8/DP/plink/varcomp.csv',
    SP8 = 'output/model8/SP/plink/varcomp.csv',
    PP8 = 'output/model8/PP/plink/varcomp.csv'
  )

  # Initialize an empty dataframe for storing summary
  summary_df <- data.frame()

  # Process each file
  for (condition in names(file_paths)) {
    # Read the data
    data <- read.csv(file_paths[[condition]], check.names = FALSE)

    # Select relevant columns
    data_selected <- data %>% select(any_of(c("V_X1", "V_X2", "V_E", "V_G", "V_GE", 'V_M7')))

    # Reshape to long format for easy summary
    data_long <- data_selected %>%
      pivot_longer(cols = everything(), names_to = "colname", values_to = "value")

    # Calculate mean and 95% CI for each component
    stats <- data_long %>%
      group_by(colname) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        ci_lower = mean(value, na.rm = TRUE) - qt(0.975, sum(!is.na(value)) - 1) * sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
        ci_upper = mean(value, na.rm = TRUE) + qt(0.975, sum(!is.na(value)) - 1) * sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
        .groups = "drop"
      ) %>%
      mutate(filename = condition) %>%
      select(filename, everything())

    # Append to final summary
    summary_df <- bind_rows(summary_df, stats)
  }

  # Save summary to file
  output_file <- "data/nklimko/plotdata/mean_ci.txt"
  write.csv(summary_df, file = output_file, row.names = FALSE)

  cat("Combined column means with confidence intervals saved to:", output_file, "\n")

  # code/scripts/nklimko/varcomp_ci.R
}




matoo <- matte[comp=='GxE',]


mathree <- matoo[order(trait, modnum),]

micols <- c('trait', 'modnum', 'mean')

mafo <- mathree[,micols,with=0]


names(mafo) <- c('Trait','Model', 'PVE')

kable(mafo, 'simple', digits=2)

names(matte)

mingus <- mat_dp[modnum > 7,]


mingus[order(comp),]



micol <- c('model','comp','mean')

tingus <- mingus[,micol,with=0]
tingus <- tingus[order(Comp),]
names(tingus) <- str_to_title(names(tingus))


kable(tingus, 'simple')

title('temp')


vignette('BGLR-extdoc', package='BGLR')




d6b <- wreath(c('DP', 'Model 6b', 'output/brit_24k/DP/plink/model6_white24k.csv'))
s6b <- wreath(c('SP', 'Model 6b', 'output/brit_24k/SP/plink/model6_white24k.csv'))
p6b <- wreath(c('PP', 'Model 6b', 'output/brit_24k/PP/plink/model6_white24k.csv'))

a6b <- rbind(d6b, s6b, p6b)

a6b$mean <- round(a6b$mean, 2)
a6b$lower <- round(a6b$lower, 2)
a6b$upper <- round(a6b$upper, 2)
a
a <- ' & '

a6b$funny <- paste0(a6b$comp, a, a6b$mean, a, a6b$lower, a, a6b$upper, ' \\ \\hline')


a6b$funny


X1 & 6.02 & 4.05 & 7.19 \\ \hline
X2 & 0.29 & 0.04 & 0.53 \\ \hline
G & 18 & 12.31 & 20.61 \\ \hline
E & 10.84 & 9.32 & 11.6 \\ \hline
GxE & 4.21 & 1.94 & 5.72 \\ \hline
X1 & 13.09 & 11.27 & 14.09 \\ \hline
X2 & 0.22 & 0.03 & 0.41 \\ \hline
G & 17.31 & 11.76 & 19.76 \\ \hline
E & 4.82 & 3.99 & 5.29 \\ \hline
GxE & 3.12 & 1.48 & 4.1 \\ \hline
X1 & 16.66 & 15.03 & 17.58 \\ \hline
X2 & 0.21 & 0.02 & 0.41 \\ \hline
G & 15.79 & 10.26 & 18.1 \\ \hline
E & 1.05 & 0.63 & 1.33 \\ \hline
GxE & 3.89 & 1.72 & 5.13 \\ \hline



