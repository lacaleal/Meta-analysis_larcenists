#### Statistical analyses - Leal et al. 
#### Submitted to Ecology

#### All packages used in the script#####

library(devtools) 
library(metafor) 
library(V.PhyloMaker) 
library(ggplot2)
library(orchaRd)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(egg)

#### Calculating the effect sizes (Hedges'g)####

data_comp=read.csv("complete_hedges.csv", sep=";", dec=)
View(data_comp)
str(data_comp)

effect_v=escalc("ROM",m1i=mean_larcenist, m2i=mean_control,
              sd1i=sd_larcenist ,sd2i=sd_control, n1i=n_larcenist,
              n2i=n_control , data=data_comp)

write.csv(effect_v, file="add_hedges.csv")

#Calling the dataset after including the columns reporting the effect sizes and its variance#

hedges=read.csv("complete_hedges.csv", sep=";", dec=".")
View(hedges)

#####Creating phylogenetic matrix - V.PhyloMaker2 package####

phylo=read.csv("lista_phylo.csv", sep=";", h=T)
View(phylo)

tree=phylo.maker(as.data.frame(phylo), tree=GBOTB.extended, nodes=nodes.info.1)
draw=write.tree(tree$scenario.3, "visit.tre")
test = read.tree("visit.tre")
a=is.ultrametric(tree$scenario.3)

#Building the covariance matrix
cov.matrix = vcv.phylo(tree$scenario.3, corr=T)
cov.matrix[,0]
View(cov.matrix)

##### Model 0: type of study #######

m0=rma.mv(yi, vi, mods= ~type_study:Larcenist_type, 
          random=list(~1|study,~1|tips_phylo,~1|sp_plant,~1|sample), 
          R = list(tips_phylo = cov.matrix), method="REML",control=list(optimizer="optim"),
          digits=3,data=hedges)
summary(m0)
confint(m0)

#Heterogeneity - I2

hO=format(i2_ml(m0, boot=100, method="ratio",
                     data=hedges))

# Publication bias #
eggsm0=lm(residuals(m0)~sqrt(as.numeric(hedges$vi))) 
summary(eggsm0)
confint(eggsm0)

# Figure S1  #

m0.1=rma.mv(yi, vi, mods= ~type_study:Larcenist_type-1, 
            random=list(~1|study,~1|tips_phylo,~1|sp_plant,~1|sample), R = list(tips_phylo = cov.matrix), method="REML",control=list(optimizer="optim"),digits=3,data=hedges)

tiff("fig_S1.tiff", units="mm", width=330, height=200, res=600, compression="lzw")

res_m0.1 <- orchaRd::mod_results(m0.1, group = "study", mod = "Larcenist_type" ,
                                 by = "type_study",
                                 data = hedges_phy2)
orchard_plot(res_m0.1,
             xlab = "Effect size (Hedges'g)",condition.lab =" Type of study",
             alpha=0.5, cb=T, k=T, trunk.size=12, branch.size = 1.5) + 
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("darkolivegreen3","orange")) +
  theme_classic() +
  theme(legend.title=element_text(size=16), legend.text= element_text(size=16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "top")
dev.off()

####Model 01 - Overall Lacernist effects####

moverall=rma.mv(yi, vi, 
                random=list(~1|study,~1|tips_phylo,~1|sp_plant,~1|sample), R = list(tips_phylo = cov.matrix), method="REML",control=list(optimizer="optim"),digits=3,data=hedges)

summary(moverall)
confint (moverall)

####Model 01 - Lacernist type ####

m1=rma.mv(yi, vi, mods= ~Larcenist_type, 
          random=list(~1|study,~1|tips_phylo,~1|sp_plant,~1|sample), 
          R = list(tips_phylo = cov.matrix), method="REML",control=list(optimizer="optim"),digits=3,data=hedges)

summary(m1)
confint(m1)

#Model heterogeneity - I2 - Orchard package

hetero1=format(i2_ml(m1, boot=100, method="ratio",
             data=hedges))

#Publication bias (egger) - 
eggsm1=lm(residuals(m1)~sqrt(as.numeric(hedges2$vi))) 
summary(eggsm1)
confint(eggsm1)

####Model 02 - Interaction between Lacernist and effect type####

m2=rma.mv(yi, vi, mods= ~Larcenist_type:Type_of_metric, 
          random=list(~1|sp_plant, ~1|tips_phylo, ~1|study,  ~1|sample), R = list(tips_phylo = cov.matrix), method="REML",control=list(optimizer="optim"),digits=3,data=hedges)
summary(m2)
confint(m2)

#Model heterogeneity - I2 - Orchard package

hetero2=format(i2_ml(m2,method="ratio",
              boot = 100, data=hedges))

#Publication bias (egger)
eggsm2=lm(residuals(m2)~sqrt(as.numeric(hedges$vi))) 
summary(eggsm2)
confint(eggsm2)

#### Data about plant performance and robbers ####

#New dataset compiling effect sizes of larcenists on plant reproductive performance

perfo_rob=read.csv("perform_robbers.csv", sep=";", h=T)
View(perfo_rob)

#### Model 03 - Matting system #####

m3=rma.mv(yi, vi, mods= ~Matting_system, 
           random=list(~1|sp_plant, ~1|tips_phylo, ~1|study,~1|sample), R = list(tips_phylo = cov.matrix), method="REML",control=list(optimizer="optim"), digits=3,data=perfo_rob)

 summary(m3)
 confint(m3)

##I-squared calculations ## 

hetero3=format(i2_ml(m3, method="ratio",
                     boot = 100, data=perfo_rob))
confint(m3)

#Publication bias (egger test)

eggsm3=lm(residuals(m3)~sqrt(as.numeric(perfo_rob$vi))) 
summary(eggsm3)
confint(eggsm3)

#### Model 04 -  Effect of players - larcenist-pollinator ID #####

m4=rma.mv(yi, vi, mods= ~players, 
          random=list( ~1|study,~1|tips_phylo,~1|sp_plant,~1|sample), 
          R = list(tips_phylo = cov.matrix), method="REML",control=list(optimizer="optim"),digits=3,data=perfo_rob)

summary(m4)
confint (m4)

#Model heterogeneity - I2 - Orchard package

hetero4=format(i2_ml(m4, boot=100, method="ratio",
                     data=perfo_rob))

#Publication bias (egger) - 

eggsm4=lm(residuals(m4)~sqrt(as.numeric(perfo_rob$vi))) 
summary(eggsm4)
confint(eggsm4)


##### Figures #####

### Fig 1 ##
m1.1=rma.mv(yi, vi, mods= ~Larcenist_type-1, 
            random=list(~1|study,~1|tips_phylo,~1|sp_plant,~1|sample), R = list(tips_phylo = cov.matrix), method="REML", digits=3, data=hedges)

res_m1 <- orchaRd::mod_results(m1.1, group = "study", mod = "Larcenist_type",
                               data = hedges)

fig1a <-orchard_plot(res_m1,mod="Larcenist_type",group="study", data=hedges_phy2,
             xlab = "Effect size (Hedges'g)", condition.lab =
               " Larcenist type", alpha=0.5, cb=T, k=T,
             trunk.size=12, branch.size = 1.5) + 
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("darkolivegreen3", "orange")) +
  theme_classic() +
  theme(legend.title=element_text(size=16), legend.text= element_text(size=16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "top")

m2.1=rma.mv(yi, vi, mods= ~Larcenist_type:Type_of_metric-1, 
            random=list(~1|sp_plant, ~1|tips_phylo, ~1|study,  ~1|sample), R = list(tips_phylo = cov.matrix), method="REML",control=list(optimizer="optim"),digits=3, data=hedges)

res_m2 <- orchaRd::mod_results(m2.1, group = "study", mod = "Type_of_metric",by = "Larcenist_type", data = hedges)

fig1b <- orchard_plot(res_m2,
                      xlab = "Effect size (Hedges'g)",condition.lab =
                        " Effect type", alpha=0.5, cb=T, k=T,
                      trunk.size=12, branch.size = 1.5)+
  scale_color_manual(values = c("black", "black", "black", "black")) +
  scale_fill_manual(values = c("darkolivegreen3","orange", "darkmagenta","cadetblue")) +
  theme_classic() +
  theme(legend.title=element_text(size=16), legend.text= element_text(size=16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "top")


tiff("fig_1_painel.tiff", units="mm", width=250, height=400, res=600, compression="lzw")

ggpubr::ggarrange(fig1a, fig1b, ncol = 1, common.legend = TRUE, labels = c("(A)", "(B)"))

dev.off ()

# Fig 2#


m3.1=rma.mv (yi, vi, mods= ~Matting_system-1, 
             random=list(~1|sp_plant, ~1|tips_phylo, ~1|study, ~1|sample),R = list(tips_phylo = cov.matrix), method="REML", digits=3,
             data=perfo_rob)

res_m3.1 <- orchaRd::mod_results(m3.1, group = "study", mod = "Matting_system",data = perfo_rob)

fig2a <- orchard_plot(res_m3.1,
                      xlab = "Effect size (Hedges'g)",condition.lab =
                        " Matting system", alpha=0.5, cb=T, k=T,
                      trunk.size=12, branch.size = 2) + 
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("darkmagenta","cadetblue")) +
  scale_x_discrete(labels = c("Selfing  ", "Not selfing  "))+
  theme_classic() +
  theme(legend.title=element_text(size=16), legend.text= element_text(size=16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "top")

# To create figure 2b, we excluded the NA values from players column in the perfo_rob dataframe and rename this new dataframe as perfo_fig

perfo_fig=read.csv("perform_rob_fig.csv", sep=";", h=T)

m4.1=rma.mv(yi, vi, mods= ~players-1, 
            random=list(~1|study,~1|tips_phylo,~1|sp_plant,~1|sample), R = list(tips_phylo = cov.matrix), method="REML",control=list(optimizer="optim"),digits=3,data=perfo_fig)

res_m4 <- orchaRd::mod_results(m4.1, group = "study", mod = "players",
                               data = perfo_fig)

fig2b <- orchard_plot(m4.1,mod="players", group="study", data=perfo_rob,xlab = "Effect size (Hedges'g)", condition.lab =
" Larcenist-pollinators", alpha=0.5, cb=T, k=T,trunk.size=10) + 
  scale_color_manual(values = c("black", "black", "black", "black", "black")) + scale_fill_manual(values = c("darkolivegreen3","darkmagenta","cadetblue")) + theme_classic() +
  theme(legend.title=element_text(size=16), legend.text= element_text(size=16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "top")

tiff("fig_2_painel.tiff", units="mm", width=250, height=400, res=600, compression="lzw")

ggpubr::ggarrange(fig2a, fig2b, ncol = 1, common.legend = TRUE, labels = c("(A)", "(B)"))

dev.off ()
