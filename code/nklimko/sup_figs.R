

library(data.table)
library(ggplot2)
library(ggpubr)


# full table of variance component effect sizes
matte <- fread('output/vault/fullport/varcomps.txt')


matter_alt <- function(datanow, fulltitle, montag=NULL, mod7=TRUE, size_point=4){
  plothole <- ggplot(data=datanow, aes(x=fax, y=mean))+
    # scale_color_manual(discrete=T, values = c('A_X1', 'B_X2', 'C_E', 'D_G', 'E_M7', 'F_GxE'), labels = c('X1', 'X2', 'E', 'G', 'M7', 'GxE'))+
    # {if(mod7)  scale_color_discrete(labels = c('X1', 'X2', 'E', 'G', 'M7', 'GxE'))}+
    # {if(!mod7) scale_color_discrete(labels = c('X1', 'X2', 'E', 'G', 'GxE'))}+
    {if(mod7)  scale_color_discrete(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Ancestry-by-lifestyle', 'Genetics-by-lifestyle'))}+
    {if(!mod7)  scale_color_discrete(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Genetics-by-lifestyle'))}+
    geom_point(aes(color=comp2), size=size_point)+
    # geom_text(aes(label=modnum))+
    geom_errorbar(aes(ymin = lower, ymax = upper, color=comp2),
                  width = 0.8, linewidth = 0.4, alpha = 1) +
    # labs(title='Model 8 Variance Components')+
    labs(title=fulltitle, color=NULL, tag=montag)+
    # xlab('Model')+
    xlab(NULL)+
    ylab('Variance')+
    theme_minimal()+
    scale_x_discrete(labels = c("Model 6_A_X1" = ' ',
                                "Model 6_B_X2" = ' ',
                                "Model 6_C_E" = 'Model 6 - Our PCs',
                                "Model 6_D_G" = ' ',
                                "Model 6_F_GxE" = ' ',
                                "Model 6a_A_X1" = ' ',
                                "Model 6a_B_X2" = ' ',
                                "Model 6a_C_E" = 'Model 6 - UKB PCs',
                                "Model 6a_D_G" = ' ',
                                "Model 6a_F_GxE" = ' '))+
    theme(legend.position = 'bottom',
          # axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    geom_vline(xintercept = verty, linetype=2, color='lightgray')+
    # geom_hline(yintercept = 0, color='black')+
    # geom_text(data=modcode, hjust='middle', aes(label=mn, x = x, y = y), inherit.aes=FALSE)+
    guides(color = guide_legend(override.aes = list(size = size_point / 2)))


  return(plothole)
}


# plotting function
matter2 <- function(datanow, fulltitle, montag=NULL, mod7=TRUE, size_point=4){
  plothole <- ggplot(data=datanow, aes(x=fax, y=mean))+
    # scale_color_manual(discrete=T, values = c('A_X1', 'B_X2', 'C_E', 'D_G', 'E_M7', 'F_GxE'), labels = c('X1', 'X2', 'E', 'G', 'M7', 'GxE'))+
    # {if(mod7)  scale_color_discrete(labels = c('X1', 'X2', 'E', 'G', 'M7', 'GxE'))}+
    # {if(!mod7) scale_color_discrete(labels = c('X1', 'X2', 'E', 'G', 'GxE'))}+
    {if(mod7)  scale_color_discrete(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Ancestry-by-lifestyle', 'Genetics-by-lifestyle'))}+
    {if(!mod7)  scale_color_discrete(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Genetics-by-lifestyle'))}+
    geom_point(aes(color=comp2), size=size_point)+
    # geom_text(aes(label=modnum))+
    geom_errorbar(aes(ymin = lower, ymax = upper, color=comp2),
                  width = 0.8, linewidth = 0.4, alpha = 1) +
    # labs(title='Model 8 Variance Components')+
    labs(title=fulltitle, color=NULL, tag=montag)+
    # xlab('Model')+
    xlab(NULL)+
    ylab('Variance')+
    theme_minimal()+
    scale_x_discrete(labels = c("Model 6_A_X1" = ' ',
                                "Model 6_B_X2" = ' ',
                                "Model 6_C_E" = 'Model 6 - Multi-ancestry',
                                "Model 6_D_G" = ' ',
                                "Model 6_F_GxE" = ' ',
                                "Model 6b_A_X1" = ' ',
                                "Model 6b_B_X2" = ' ',
                                "Model 6b_C_E" = 'Model 6 - British',
                                "Model 6b_D_G" = ' ',
                                "Model 6b_F_GxE" = ' '))+
    theme(legend.position = 'bottom',
          # axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    geom_vline(xintercept = verty, linetype=2, color='lightgray')+
    # geom_hline(yintercept = 0, color='black')+
    # geom_text(data=modcode, hjust='middle', aes(label=mn, x = x, y = y), inherit.aes=FALSE)+
    guides(color = guide_legend(override.aes = list(size = size_point / 2)))


  return(plothole)
}


# plotting function
matter7 <- function(datanow, fulltitle, montag=NULL, mod7=TRUE, size_point=4, mjust=0, textjust=3){
  plothole <- ggplot(data=datanow, aes(x=fax, y=mean))+
    {if(mod7)  scale_color_discrete(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Ancestry-by-lifestyle', 'Genetics-by-lifestyle'))}+
    {if(!mod7)  scale_color_discrete(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Genetics-by-lifestyle'))}+
    geom_point(aes(color=comp2), size=size_point, )+
    geom_text(aes(label =round(mean, 2), x=fax, y=mean + textjust), inherit.aes = FALSE)+
    geom_errorbar(aes(ymin = lower, ymax = upper, color=comp2),
                  width = 0.8, linewidth = 0.4, alpha = 1) +
    labs(title=fulltitle, color=NULL, tag=montag)+
    xlab(NULL)+
    ylab('Variance')+
    theme_minimal()+
    scale_x_discrete(labels = c("Model 7_A_X1" = ' ',
                                "Model 7_B_X2" = ' ',
                                "Model 7_C_E" = 'Model 7 - Ancestry-by-lifestyle',
                                "Model 7_D_G" = ' ',
                                "Model 7_E_M7" = ' ',
                                "Model 7_F_GxE" = ' '))+
    theme(legend.position = 'bottom',
          axis.text.x = element_text(hjust = mjust),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    geom_vline(xintercept = verty, linetype=2, color='lightgray')+
    # geom_hline(yintercept = 0, color='black')+
    # geom_text(data=modcode, hjust='middle', aes(label=mn, x = x, y = y), inherit.aes=FALSE)+
    guides(color = guide_legend(override.aes = list(size = size_point / 2)))


  return(plothole)
}




plot_d7 <- matter7(mat_dp7, fulltitle=NULL, montag='A', mod7=TRUE, mjust=0.1)
plot_s7 <- matter7(mat_sp7, fulltitle=NULL, montag='B', mod7=TRUE, mjust=0.1)
plot_p7 <- matter7(mat_pp7, fulltitle=NULL, montag='C', mod7=TRUE, mjust=0.1)


triple <- ggarrange(plot_d7, plot_s7, plot_p7, ncol=1, common.legend = TRUE, legend='bottom')

triple


matter7(mat_dp7, fulltitle=NULL, montag='A', mod7=TRUE, mjust=0.1, textjust=1)
#
# plot_d7


matterall <- function(datanow, fulltitle, montag=NULL, mod7=TRUE, size_point=4, mjust=0){
  plothole <- ggplot(data=datanow, aes(x=fax, y=mean))+
    {if(mod7)  scale_color_discrete(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Ancestry-by-lifestyle', 'Genetics-by-lifestyle'))}+
    {if(!mod7)  scale_color_discrete(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Genetics-by-lifestyle'))}+
    geom_point(aes(color=comp2), size=size_point)+
    geom_errorbar(aes(ymin = lower, ymax = upper, color=comp2),
                  width = 0.8, linewidth = 0.4, alpha = 1) +
    labs(title=fulltitle, color=NULL, tag=montag)+
    xlab(NULL)+
    ylab('Variance')+
    theme_minimal()+
    scale_x_discrete(labels = c("Model 1_A_X1" = 'Model 1',
                                "Model 2_A_X1" = 'Model 2',
                                "Model 2_B_X2" = ' ',
                                "Model 3_A_X1" = 'Model 3',
                                "Model 3_B_X2" = ' ',
                                "Model 3_D_G" = ' ',
                                "Model 4_A_X1" = 'Model 4',
                                "Model 4_B_X2" = ' ',
                                "Model 4_C_E" = ' ',
                                "Model 5_A_X1" = 'Model 5',
                                "Model 5_B_X2" = ' ',
                                "Model 5_C_E" = ' ',
                                "Model 5_D_G" = ' ',
                                "Model 6_A_X1" = 'Model 6',
                                "Model 6_B_X2" = ' ',
                                "Model 6_C_E" = ' ',
                                "Model 6_D_G" = ' ',
                                "Model 6_F_GxE" = ' ',
                                "Model 7_A_X1" = 'Model 7',
                                "Model 7_B_X2" = ' ',
                                "Model 7_C_E" = ' ',
                                "Model 7_D_G" = ' ',
                                "Model 7_E_M7" = ' ',
                                "Model 7_F_GxE" = ' '))+
    theme(legend.position = 'bottom',
          # axis.text.x = element_text(hjust = mjust),
          # axis.text.x = element_text(size=8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    geom_vline(xintercept = verty, linetype=2, color='lightgray')+
    # geom_hline(yintercept = 0, color='black')+
    # geom_text(data=modcode, hjust='middle', aes(label=mn, x = x, y = y), inherit.aes=FALSE)+
    guides(color = guide_legend(override.aes = list(size = size_point / 2)))


  return(plothole)
}

tempy <- matterall(mat_dp_num, fulltitle=NULL, montag='DP', mod7=TRUE, mjust = NA)
triple <- ggarrange(tempy, tempy, tempy, ncol=1, common.legend = TRUE, legend='bottom')

triple

# scale_x_discrete(labels = c())

scale_x_discrete(labels = c("Model 1_A_X1" = 'Model 1',
                            "Model 2_A_X1" = 'Model 2',
                            "Model 2_B_X2" = ' ',
                            "Model 3_A_X1" = 'Model 3',
                            "Model 3_B_X2" = ' ',
                            "Model 3_D_G" = ' ',
                            "Model 4_A_X1" = 'Model 4',
                            "Model 4_B_X2" = ' ',
                            "Model 4_C_E" = ' ',
                            "Model 5_A_X1" = 'Model 5',
                            "Model 5_B_X2" = ' ',
                            "Model 5_C_E" = ' ',
                            "Model 5_D_G" = ' ',
                            "Model 6_A_X1" = 'Model 6',
                            "Model 6_B_X2" = ' ',
                            "Model 6_C_E" = ' ',
                            "Model 6_D_G" = ' ',
                            "Model 6_F_GxE" = ' ',
                            "Model 6a_A_X1" = 'Model 6 - UK Biobank PCs',
                            "Model 6a_B_X2" = ' ',
                            "Model 6a_C_E" = ' ',
                            "Model 6a_D_G" = ' ',
                            "Model 6a_F_GxE" = ' ',
                            "Model 7_A_X1" = 'Model 7',
                            "Model 7_B_X2" = ' ',
                            "Model 7_C_E" = ' ',
                            "Model 7_D_G" = ' ',
                            "Model 7_E_M7" = ' ',
                            "Model 7_F_GxE" = ' '))

# Full read in
if(1){


  # Variance component readin ----

  # split by trait
  mat_dp <- matte[trait=='DP',]
  mat_sp <- matte[trait=='SP',]
  mat_pp <- matte[trait=='PP',]

  # alt model 6 gen pcs
  mod1 <- c('6', '6a')
  mat_dp_pcs <- mat_dp[modnum %in% mod1,]
  mat_sp_pcs <- mat_sp[modnum %in% mod1,]
  mat_pp_pcs <- mat_pp[modnum %in% mod1,]

  # alt model 6 ancestry grouping
  mod2 <- c('6', '6b')
  # mod2 <- c('6', '6b', '6w')
  mat_dp_alt <- mat_dp[modnum %in% mod2,]
  mat_sp_alt <- mat_sp[modnum %in% mod2,]
  mat_pp_alt <- mat_pp[modnum %in% mod2,]

  # base model 6 vs new models 7-9
  mod3 <- as.character(1:7)
  mat_dp_num <- mat_dp[modnum %in% mod3,]
  mat_sp_num <- mat_sp[modnum %in% mod3,]
  mat_pp_num <- mat_pp[modnum %in% mod3,]

  # model 7 alone
  mod4 <- as.character(7)
  mat_dp7 <- mat_dp[modnum %in% mod4,]
  mat_sp7 <- mat_sp[modnum %in% mod4,]
  mat_pp7 <- mat_pp[modnum %in% mod4,]

  # Variance component plot formatting ----
  # indexes of line breaks - model 6 pcs
  # verty <- c(1.5, 3.5, 6.5, 9.5, 13.5, 17.5, 22.5, 27.5, 33.5)
  verty <- 0.5 + cumsum(unlist(table(mat_dp_pcs$model)))
  verty <- verty[-length(verty)]

  # add on softpoint ends
  verty2 <- c(0.5, verty, nrow(mat_dp_pcs)+0.5)

  # calculates midpoint of each section to place Model Number below x axis
  verty3 <- (verty2[2:length(verty2)] - verty2[1:(length(verty2) - 1)]) / 2 + verty2[1:(length(verty2) - 1)]

  altnumber <- c('6 - Original', '6a - Alt PCs')

  #matter function plots model number at specified breakpoint at y = -1 for each label.
  modcode <- data.table(mn=altnumber, x=verty3, y=-1)

  dpc <- matter_alt(mat_dp_pcs, fulltitle=NULL, montag='DP', mod7=FALSE)
  spc <- matter_alt(mat_sp_pcs, fulltitle=NULL, montag='SP', mod7=FALSE)
  ppc <- matter_alt(mat_pp_pcs, fulltitle=NULL, montag='PP', mod7=FALSE)

  print(plot_grid(dpc, spc, ppc, ncol=1))

  # indexes of line breaks - model 6 white groups
  # verty <- c(1.5, 3.5, 6.5, 9.5, 13.5, 17.5, 22.5, 27.5, 33.5)

  mat_dp_alt <- mat_dp_alt[modnum!='6w',]
  mat_sp_alt <- mat_sp_alt[modnum!='6w',]
  mat_pp_alt <- mat_pp_alt[modnum!='6w',]

  verty <- 0.5 + cumsum(unlist(table(mat_dp_alt$model)))
  verty <- verty[-length(verty)]

  dalt <- matter2(mat_dp_alt, fulltitle=NULL, montag='A', mod7=FALSE, size_point=4)
  salt <- matter2(mat_sp_alt, fulltitle=NULL, montag='B', mod7=FALSE, size_point=4)
  palt <- matter2(mat_pp_alt, fulltitle=NULL, montag='C', mod7=FALSE, size_point=4)

  # print(plot_grid(dalt, salt, palt, ncol=1))

  # indexes of line breaks - models 7-9
  # verty <- c(1.5, 3.5, 6.5, 9.5, 13.5, 17.5, 22.5, 27.5, 33.5)
  verty <- 0.5 + cumsum(unlist(table(mat_dp_num$model)))
  verty <- verty[-length(verty)]

  dnum <- matterall(mat_dp_num, fulltitle=NULL, montag='DP', mod7=TRUE)
  snum <- matterall(mat_sp_num, fulltitle=NULL, montag='SP', mod7=TRUE)
  pnum <- matterall(mat_pp_num, fulltitle=NULL, montag='PP', mod7=TRUE)

  verty <- 0.5 + cumsum(unlist(table(mat_dp7$model)))
  verty <- verty[-length(verty)]

  plot_d7 <- matter7(mat_dp7, fulltitle=NULL, montag='A', mod7=TRUE, mjust=0.1)
  plot_s7 <- matter7(mat_sp7, fulltitle=NULL, montag='B', mod7=TRUE, mjust=0.1)
  plot_p7 <- matter7(mat_pp7, fulltitle=NULL, montag='C', mod7=TRUE, mjust=0.1)



  verty <- 0.5 + cumsum(unlist(table(mat_dp_num$model)))
  verty <- verty[-length(verty)]



  dnum2 <- matterall2(mat_dp_num, fulltitle=NULL, montag='A', mod7=TRUE)
  snum2 <- matterall(mat_sp_num, fulltitle=NULL, montag='B', mod7=TRUE)
  pnum2 <- matterall(mat_pp_num, fulltitle=NULL, montag='C', mod7=TRUE)



}
#
#   ggplot(df %>% filter(trait == trait_name),
#          aes(x = Model_simple, y = mean, color = Component_label)) +
#     geom_point(position = position_dodge2(width = 1, preserve = "single"), size = 9) +
#     geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 1,
#                   position = position_dodge2(width = 1, preserve = "single")) +

matterall2 <- function(datanow, fulltitle, montag=NULL, mod7=TRUE, size_point=4, mjust=0, textjust=0, buffer=0, custard = 'nerg'){
  plothole <- ggplot(data=datanow, aes(x=modflex, y=mean))+
    {if(mod7)  scale_color_manual(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Ancestry-by-lifestyle', 'Genetics-by-lifestyle', ' '), values= c("#E69F00","#56B4E9","#009E73","purple",'pink', "brown",'black'))}+
    {if(!mod7)  scale_color_manual(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Genetics-by-lifestyle', ' '), values= c("#E69F00","#56B4E9","#009E73","purple","brown",'white'))}+
    # geom_point(aes(color=comp2), size=size_point)+
    geom_point(aes(color=comp2), position = position_dodge2(width = 1, preserve = "single"), size = size_point/2) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color=comp2),
                  width = 1, linewidth = 0.4, alpha = 1, position = position_dodge2(width = 1, preserve = "single")) +
    # geom_text(aes(label =round(mean, 2), y = textjust, color=comp2), angle=90,  size=size_point - 1 ,position = position_dodge2(width = 1, preserve = "single"))+
    # geom_text(aes(label = custard, y = textjust, color='black'), size=size_point , position = position_dodge2(width = 1, preserve = "single"))+
    # annotate(geom_text, x = 0, y = textjust + 2, color='white', label='.')
    geom_text(data=custard, aes(x = modflex, y = textjust, label= round(V1, 2)), color='black', show.legend = FALSE)+
    geom_text(aes(x = 'Model 2', y = textjust + buffer, label=" ", color='white', fill='white'), show.legend = FALSE)+
    labs(title=fulltitle, color=NULL, tag=montag)+
    xlab(NULL)+
    ylab('Variance')+
    theme_minimal()+
    scale_x_discrete(labels = c("Model 1_A_X1" = 'Model 1',
                                "Model 2_A_X1" = 'Model 2',
                                "Model 2_B_X2" = ' ',
                                "Model 3_A_X1" = 'Model 3',
                                "Model 3_B_X2" = ' ',
                                "Model 3_D_G" = ' ',
                                "Model 4_A_X1" = 'Model 4',
                                "Model 4_B_X2" = ' ',
                                "Model 4_C_E" = ' ',
                                "Model 5_A_X1" = 'Model 5',
                                "Model 5_B_X2" = ' ',
                                "Model 5_C_E" = ' ',
                                "Model 5_D_G" = ' ',
                                "Model 6_A_X1" = 'Model 6',
                                "Model 6_B_X2" = ' ',
                                "Model 6_C_E" = ' ',
                                "Model 6_D_G" = ' ',
                                "Model 6_F_GxE" = ' ',
                                "Model 7_A_X1" = 'Model 7',
                                "Model 7_B_X2" = ' ',
                                "Model 7_C_E" = ' ',
                                "Model 7_D_G" = ' ',
                                "Model 7_E_M7" = ' ',
                                "Model 7_F_GxE" = ' '))+
    theme(legend.position = 'bottom',
          # axis.text.x = element_text(hjust = mjust),
          # axis.text.x = element_text(size=8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    # geom_vline(xintercept = verty, linetype=2, color='lightgray')+
    # geom_hline(yintercept = 0, color='black')+
    # geom_text(data=modcode, hjust='middle', aes(label=mn, x = x, y = y), inherit.aes=FALSE)+
    guides(color = guide_legend(override.aes = list(size = size_point / 2)))


  return(plothole)
}

matterall2(mat_dp_num, fulltitle=NULL, montag='A', mod7=TRUE, textjust=21, buffer = 3, custard=dp_bind)

mat_dp_num$modflex <- paste0('Model ',mat_dp_num$modnum)
mat_sp_num$modflex <- paste0('Model ',mat_sp_num$modnum)
mat_pp_num$modflex <- paste0('Model ',mat_pp_num$modnum)

keybind <- data.table(paste0('Model ',1:7), 1:7)


dp_bind <- mat_dp_num[,sum(mean), by=.(modflex)]
sp_bind <- mat_sp_num[,sum(mean), by=.(modflex)]
pp_bind <- mat_pp_num[,sum(mean), by=.(modflex)]

dpom <- matterall2(mat_dp_num, fulltitle=NULL, montag='A', mod7=TRUE, textjust=21, buffer=2, custard=dp_bind)
spom <- matterall2(mat_sp_num, fulltitle=NULL, montag='B', mod7=TRUE, textjust=21, buffer=2, custard=sp_bind)
ppom <- matterall2(mat_pp_num, fulltitle=NULL, montag='C', mod7=TRUE, textjust=25, buffer=2, custard=pp_bind)

if(1){

  triple <- ggarrange(plot_d7, plot_s7, plot_p7, ncol=1, common.legend = TRUE, legend='bottom')

  root <- 'output/vault/sups/'
  myname <- 'm7'
  ext <- 'pdf'
  dpi <- 600

  munit <- 'in'
  finH <- 7
  finW <- 8

  myfile <- paste0(myname, '.', ext)
  fullfile <- paste0(root, myfile)

  ggsave(triple, path=root, filename=myfile, device=ext, units=munit, height = finH, width=finW)

}

# Full save sups_alt6.pdf
if(1){
  # Save final figure

  triple <- ggarrange(plot_d7, plot_s7, plot_p7, ncol=1, common.legend = TRUE, legend='bottom')

  root <- 'output/vault/sups/'
  myname <- 'm7'
  ext <- 'pdf'
  dpi <- 600

  munit <- 'in'
  finH <- 7
  finW <- 8

  myfile <- paste0(myname, '.', ext)
  fullfile <- paste0(root, myfile)

  ggsave(triple, path=root, filename=myfile, device=ext, units=munit, height = finH, width=finW)


}

# Model 7 alone
# Full save m7.pdf
if(1){

  #   plot_d7 <- matter7(mat_dp7, fulltitle=NULL, montag='A', mod7=TRUE, mjust=0.1)
  # plot_s7 <- matter7(mat_sp7, fulltitle=NULL, montag='B', mod7=TRUE, mjust=0.1)
  # plot_p7 <- matter7(mat_pp7, fulltitle=NULL, montag='C', mod7=TRUE, mjust=0.1)


  triple <- ggarrange(plot_d7, plot_s7, plot_p7, ncol=1, common.legend = TRUE, legend='bottom')

  root <- 'output/vault/sups/'
  myname <- 'm7'
  ext <- 'pdf'
  dpi <- 600

  munit <- 'in'
  finH <- 7
  finW <- 8

  myfile <- paste0(myname, '.', ext)
  fullfile <- paste0(root, myfile)

  ggsave(triple, path=root, filename=myfile, device=ext, units=munit, height = finH, width=finW)

}

# Models 1 through 7
# Full save all_models.pdf
if(1){

  # triple <- ggarrange(dnum, snum , pnum, ncol=1, common.legend = TRUE, legend='bottom')
  triple <- ggarrange(dpom, spom, ppom, ncol=1, common.legend = TRUE, legend='bottom')

  root <- 'output/vault/sups/'
  # myname <- 'all_models'
  myname <- 'fig4_final'
  ext <- 'pdf'
  dpi <- 600

  munit <- 'in'
  finH <- 7
  finW <- 8

  myfile <- paste0(myname, '.', ext)
  fullfile <- paste0(root, myfile)

  ggsave(triple, path=root, filename=myfile, device=ext, units=munit, height = finH, width=finW)

}



# Color palette for components (Genetics and Lifestyle colors swapped)
color_blind_palette <- c(
  "Demographics" = "#E69F00",
  "Structure" = "#56B4E9",
  "Genetics" = "purple",  # Swapped color
  "Lifestyle" = "#009E73", # Swapped color
  "Genetics-by-lifestyle" = "brown",
  "Ancestry-by-lifestyle" = 'pink'
)







matter <- function(datanow, fulltitle, montag=NULL){
  plothole <- ggplot(data=datanow, aes(x=fax, y=mean))+
    # scale_color_manual(discrete=T, values = c('A_X1', 'B_X2', 'C_E', 'D_G', 'E_M7', 'F_GxE'), labels = c('X1', 'X2', 'E', 'G', 'M7', 'GxE'))+
    scale_color_discrete(labels = c('X1', 'X2', 'E', 'G', 'M7', 'GxE'))+
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








matterall3 <- function(datanow, fulltitle, montag=NULL, mod7=TRUE, size_point=4, mjust=0, textjust=0, buffer=0, custard = 'nerg'){
  plothole <- ggplot(data=datanow, aes(x=modflex, y=mean))+
    {if(mod7)  scale_color_manual(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Ancestry-by-lifestyle', 'Genetics-by-lifestyle', ' '), values= c("#E69F00","#56B4E9","#009E73","purple",'pink', "brown",'black'))}+
    {if(!mod7)  scale_color_manual(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Genetics-by-lifestyle', ' '), values= c("#E69F00","#56B4E9","#009E73","purple", "brown",'black'))}+
    # geom_point(aes(color=comp2), size=size_point)+
    geom_point(aes(color=comp2), position = position_dodge2(width = 1, preserve = "single"), size = size_point/2) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color=comp2),
                  width = 1, linewidth = 0.4, alpha = 1, position = position_dodge2(width = 1, preserve = "single")) +
    # geom_text(aes(label =round(mean, 2), y = textjust, color=comp2), angle=90,  size=size_point - 1 ,position = position_dodge2(width = 1, preserve = "single"))+
    # geom_text(aes(label = custard, y = textjust, color='black'), size=size_point , position = position_dodge2(width = 1, preserve = "single"))+
    # annotate(geom_text, x = 0, y = textjust + 2, color='white', label='.')
    geom_text(data=custard, aes(x = modflex, y = textjust, label= round(V1, 2)), color='black', show.legend = FALSE)+
    geom_text(aes(x = 0, y = textjust + buffer, label=".", color='white'), show.legend = FALSE)+
    labs(title=fulltitle, color=NULL, tag=montag)+
    xlab(NULL)+
    ylab('Variance')+
    theme_minimal()+
    scale_x_discrete(labels = c("Model 6_A_X1" = 'Model 6',
                                "Model 6_B_X2" = ' ',
                                "Model 6_C_E" = ' ',
                                "Model 6_D_G" = ' ',
                                "Model 6_F_GxE" = ' ',
                                "Model 6a_A_X1" = 'Model 6 alt',
                                "Model 6a_B_X2" = ' ',
                                "Model 6a_C_E" = ' ',
                                "Model 6a_D_G" = ' ',
                                "Model 6a_F_GxE" = ' '))+
    theme(legend.position = 'bottom',
          # axis.text.x = element_text(hjust = mjust),
          # axis.text.x = element_text(size=8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    # geom_vline(xintercept = verty, linetype=2, color='lightgray')+
    # geom_hline(yintercept = 0, color='black')+
    # geom_text(data=modcode, hjust='middle', aes(label=mn, x = x, y = y), inherit.aes=FALSE)+
    guides(color = guide_legend(override.aes = list(size = size_point / 2)))


  return(plothole)
}



if(sub1){
  sub1 <- 1
  # indexes of line breaks
  # verty <- c(1.5, 3.5, 6.5, 9.5, 13.5, 17.5, 22.5, 27.5, 33.5)
  verty <- 0.5 + cumsum(unlist(table(mat_dp_pcs$model)))
  verty <- verty[-length(verty)]

  # add on hardpoint ends
  verty2 <- c(0, verty, nrow(mat_dp_pcs)+1)

  # calculates midpoint of each section to place Model Number below x axis
  verty3 <- (verty2[2:length(verty2)] - verty2[1:(length(verty2) - 1)]) / 2 + verty2[1:(length(verty2) - 1)]

  altnumber <- c('6 - Original', '6a - Alt PCs')

  #matter function plots model number at specified breakpoint at y = -1 for each label.
  modcode <- data.table(mn=altnumber, x=verty3, y=-1)

  dpc <- matterall3(mat_dp_pcs, fulltitle=NULL, montag='DP', mod7=FALSE)
  spc <- matterall3(mat_sp_pcs, fulltitle=NULL, montag='SP', mod7=FALSE)
  ppc <- matterall3(mat_pp_pcs, fulltitle=NULL, montag='PP', mod7=FALSE)



mat_dp_pcs$modflex <- paste0('Model ',mat_dp_pcs$modnum)
mat_sp_pcs$modflex <- paste0('Model ',mat_sp_pcs$modnum)
mat_pp_pcs$modflex <- paste0('Model ',mat_pp_pcs$modnum)

keybind <- data.table(paste0('Model ',1:7), 1:7)


dp_bind <- mat_dp_pcs[,sum(mean), by=.(modflex)]
sp_bind <- mat_sp_pcs[,sum(mean), by=.(modflex)]
pp_bind <- mat_pp_pcs[,sum(mean), by=.(modflex)]

dpc <- matterall2(mat_dp_pcs, fulltitle=NULL, montag='A', mod7=FALSE, textjust=21, buffer=2, custard=dp_bind)
spc <- matterall2(mat_sp_pcs, fulltitle=NULL, montag='B', mod7=FALSE, textjust=21, buffer=2, custard=sp_bind)
ppc <- matterall2(mat_pp_pcs, fulltitle=NULL, montag='C', mod7=FALSE, textjust=25, buffer=2, custard=pp_bind)


  print(plot_grid(dpc, spc, ppc, ncol=1))

}



if(1){

  triple <- ggarrange(dpc, spc, ppc, ncol=1, common.legend = TRUE, legend='bottom')

  root <- 'output/vault/sups/'
  myname <- 'm6_pcs_new'
  ext <- 'pdf'
  dpi <- 600

  munit <- 'in'
  finH <- 7
  finW <- 8

  myfile <- paste0(myname, '.', ext)
  fullfile <- paste0(root, myfile)

  ggsave(triple, path=root, filename=myfile, device=ext, units=munit, height = finH, width=finW)

}
