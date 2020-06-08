##########################################################################
### Wing Shape Integration Measures, wing defects and Figure 3, Supplementary Figures 2 and 6 script ####
##########################################################################

########## Load Libraries and local functions #################
library(effects)
library(ggplot2)
library(lme4)
library(emmeans)
library(glmmTMB)
library(tidyr)
library(dplyr)
library(geomorph)
library(car)
library(reshape2)
library(MCMCglmm)
library(RRPP)
library(boot)
library(cowplot)
library(magick)
library(gtable)
library(grid)
library(gridExtra)
library(multipanelfigure)
library(egg)
library(patchwork)
###############################################################

######################## Load Data ############################

source('../scripts/MP_African_ShortDataCleaning.R')

###############################################################


#################### Set ggplot theme #########################

pref_theme <- theme_classic() + 
  theme(text = element_text(size = 65), 
        axis.title.y = element_text(margin = margin(t = 10, r = 20, b = 10, l = 10)),
        axis.title.x = element_text(margin = margin(t = 20, r = 10, b = 10, l = 10)),
        axis.text = element_text(margin = margin(t = 15, r = 15, b = 15, l = 15)))
theme_set(pref_theme)

###############################################################

############# Calculating Wing Shape Tot Var Eccentricty rSDE and rSDE ##############
#PCA
wing_PCs <- prcomp(wings[,10:105])
wingsPC <- data.frame(wings, wing_PCs$x[,1:58])

# From Annat's code (Sept 2016)
# Getting rid of sex variation? (IAN AM I RIGHT?)
Xall.res <- matrix(NA, nrow = nrow(wingsPC), ncol = 58) # ncol is the number of PCs we want to use
line_values <- rep(NA, length = nrow(wingsPC)) #  A vector with line names.

for (ln in levels(wingsPC$line)) {
  ind <- which(wingsPC$line==ln)
  Xln <- as.matrix(wingsPC[ind,110:167]) # extracting out the right rows (line) and columns
  lcs <- log(wingsPC$CS_rescaled[ind]) # log centroid size
  sex <- wingsPC$sex[ind]
  reg <- lm(Xln ~ sex + lcs) # regression model
  Xm <- colMeans(predict(reg)) #CHECK IF ITS PREDICTING THE RIGHT WAY # predicted mean configuration at mean centroid size to add back to the residuals. Predict is fitting it by default at the mean for lcs.
  Xall.res[ind,] <- reg$residuals + matrix(Xm, nr=length(ind), nc=length(Xm), byrow=TRUE) # recentered on the line mean
  line_values[ind] = ln
}

# line_values <- as.factor(line_values) # just a double check that the order is maintained.
# data.frame(line_values, wingsPC$line)

# Making a dataframe of centered shape
centered.shape <- data.frame(wingsPC$line, Xall.res)

# Producing empty vectors for total variance, eccentricity rSDE and rSDE2 calculations
tot_var <- rep(NA, nlevels(centered.shape$wingsPC.line))
eccentricity <- rep(NA, nlevels(centered.shape$wingsPC.line))
rSDE <- rep(NA, nlevels(centered.shape$wingsPC.line))
rSDE2 <- rep(NA, nlevels(centered.shape$wingsPC.line))


# Calculating total variance, eccentricity rSDE and rSDE2 values
for (ln in 1:nlevels(centered.shape$wingsPC.line)) {
  lev <- levels(centered.shape$wingsPC.line)[ln]
  mat = centered.shape[centered.shape$wingsPC.line == lev, 2:59]
  p <- ncol(mat) # number of variables
  cov_mat <- cov(mat)
  eig_mat <- eigen(cov_mat)$values
  tot_var[ln] <- sum(eig_mat) 
  
  eccentricity[ln] <- eig_mat[1]/(sum(eig_mat))
  
  rSDE[ln] <- sqrt(var(eig_mat)*((p-1)/p)) # rescaled to p instead of p - 1 
  
  
  rSDE2[ln] <- sqrt(var(eig_mat/sum(eig_mat))*((p-1)/p) )# rescaling this by the total variance. I think since we have already applied some corrections for size and sex this may result in over-correction.
}

#making a data frame of total variance, eccentricity rSDE and rSDE2 values, and clean up
integration_measures <- data.frame(tot_var, eccentricity, rSDE, rSDE2, levels(centered.shape$wingsPC.line))

integration_measures$population <- as.factor(c(rep("Ethiopia", 27), rep("Zambia", 20)))
#############################################################

############ Integration Measures Analysis ##################

## Total Variance GLMM model

# multiplying total variance values by 1000 to make them more plot friendly and easier for model to handle
integration_measures$tot_var_alt <- integration_measures$tot_var * 1000

names(integration_measures) <- c("tot_var", "eccentricity", "rSDE", "rSDE2", 
                                 "line", "population", "tot_var_alt")

lmm.TotVar.gamma <- glmmTMB(tot_var_alt ~ population + (1|line:population),
                                 data = integration_measures, family = "Gamma")

summary(lmm.TotVar.gamma)

#Table S9: Total Variance ANOVA
Anova(lmm.TotVar.gamma)

# Extracting effects for Total variance model

lmm.TotVar.gamma_marginal_means <- as.data.frame(emmeans(lmm.TotVar.gamma, ~ population ))

#Figure 3A: Total Variance by population
n <- ggplot(lmm.TotVar.gamma_marginal_means, aes( y = (1/emmean), x = population, 
                                                        colour=population, shape = population)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=(1/lower.CL), ymax= (1/upper.CL)), 
                lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = integration_measures, aes(x = population, y = tot_var_alt),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Wing Shape Variation (TV)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
n


## Eccentricity GLMM model
lmm.Ecc.gamma <- glmmTMB(eccentricity ~ population,
                        data = integration_measures, family = "Gamma")

summary(lmm.Ecc.gamma)

# Table S10: Eccentricity ANOVA
Anova(lmm.Ecc.gamma)

# Extracting effects from Eccentricity Model
lmm.Ecc.gamma_marginal_means <- as.data.frame(emmeans(lmm.Ecc.gamma, ~ population ))


# Figure 3B: Eccentricity by population
o <- ggplot(lmm.Ecc.gamma_marginal_means, aes( y = (1/emmean), x = population, 
                                                    colour=population, shape = population)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=(1/lower.CL), ymax= (1/upper.CL)), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = integration_measures, aes(x = population, y = eccentricity),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Wing Shape Intergration (Ecc)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
o

## rSDE GLMM model

# Multiplying rSDE by 10000 to make the values easier to plot and model
integration_measures$rSDE_alt <- integration_measures$rSDE * 10000


lmm.rSDE.gamma <- glmmTMB(rSDE_alt ~ population,
                     data = integration_measures, family = "Gamma")

summary(lmm.rSDE.gamma)

# Table S11: rSDE ANOVA
Anova(lmm.rSDE.gamma)

#Extracting rSDE model effects
lmm.rSDE.gamma_marginal_means <- as.data.frame(emmeans(lmm.rSDE.gamma, ~ population ))


p <- ggplot(lmm.rSDE.gamma_marginal_means, aes( y = (1/emmean), x = population, 
                                                 colour=population, shape = population)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=(1/lower.CL), ymax=(1/upper.CL) ), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = integration_measures, aes(x = population, y = rSDE_alt),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Wing Shape Integration (rSDE)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
p


## rSDE2 GLMM model

lmm.rSDE2.gamma <- glmmTMB(rSDE2 ~ population,
                      data = integration_measures, family = "Gamma")

summary(lmm.rSDE2.gamma)

# Table 12: rSDE2 ANOVA
Anova(lmm.rSDE2.gamma)

#Extracting rSDE2 model effects
lmm.rSDE2.gamma_marginal_means <- as.data.frame(emmeans(lmm.rSDE2.gamma, ~ population ))


q <- ggplot(lmm.rSDE2.gamma_marginal_means, aes( y = (1/emmean), x = population, 
                                                  colour=population, shape = population)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=(1/lower.CL), ymax= (1/upper.CL) ), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = integration_measures, aes(x = population, y = rSDE2),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Wing Shape Integration (rSDE2)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
q


# for legends
q_2 <- ggplot(integration_measures, aes( y = rSDE2, x = population, color = population, shape = population)) +
  geom_boxplot(show.legend= FALSE, lwd=1, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), size = 10, show.legend= TRUE)+
  theme(legend.position = c(0.5, 0.9)) +
  labs(x="", y="Wing Shape Variation (rSDE2)")+
  scale_x_discrete(name="", 
                   labels = c("High-Altitude","Low-Altitude"))+
  scale_color_manual(name = "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60")) +
  scale_shape_discrete(name = "Population", 
                       labels = c("High-Altitude","Low-Altitude"))

q_2
###############################################################

############# Wing Shape Tot Var and Ecc vs mutants ##############

# Loading and cleaning mutant data
mutant_frequencies <-read.csv("../data/African_Mutation_frequencies.csv", h = T)

mutant_frequencies$alv_prop <- with(mutant_frequencies, alv_num/total_wings)
mutant_frequencies$micv_prop <- with(mutant_frequencies, micv_num/total_wings)
mutant_frequencies$defects_prop <- with(mutant_frequencies, total_defects/total_wings)

# create a population variable 
boolean_for_line <- grepl(pattern = "EF", mutant_frequencies$line)
mutant_frequencies$population <- as.factor(ifelse(boolean_for_line, "Ethiopia", "Zambia"))

integration_measures2 <- integration_measures 


integration_measures2$line <- as.factor(toupper(integration_measures2$line)) 
# average mutant phenotypes across sex
mutant_frequencies_averaged <- aggregate(x = mutant_frequencies[,c(11, 12, 13 )], 
                                         by = c(mutant_frequencies["population"], mutant_frequencies["line"]), 
                                         FUN = mean )

# Merging defects and integration measures data
mutant_integration <- merge(x = integration_measures2, y = mutant_frequencies_averaged, by = c("line", "population") ) 

# Supplementary Figure S2A: Total Variance vs Proportion of defects
r <- ggplot(mutant_integration, aes( y =defects_prop, x = tot_var_alt , 
                                     color = population, shape = population)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.8, 0.9), 
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  scale_x_continuous(name="Wing Shape Variation (TV)")+
  scale_y_continuous(name="Proportion of Defects")+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population",
                      labels = c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))
r 

# Correlation tests for Total variance vs wing defects
mutant_integration_Eth <- subset(mutant_integration, population == "Ethiopia")

cor.test(mutant_integration_Eth$tot_var_alt, mutant_integration_Eth$defects_prop)

mutant_integration_Zam <- subset(mutant_integration, population == "Zambia")

cor.test(mutant_integration_Zam$tot_var_alt, mutant_integration_Zam$defects_prop)


# Supplementary Figure S2B: Eccentricity vs Proportion of defects
s <- ggplot(mutant_integration, aes( y = defects_prop, x =eccentricity , 
                                     color = population, shape = population)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.8, 0.9), 
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  scale_x_continuous(name="Wing Shape Integration (Ecc)")+
  scale_y_continuous(name="Proportion of Defects")+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population",
                      labels = c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))
s

s_2 <- ggplot(mutant_integration, aes( y = defects_prop, x =eccentricity , 
                                     color = population, shape = population)) +
  geom_point(size = 10, show.legend = TRUE)+
  theme(legend.position = c(0.8, 0.9), 
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  scale_x_continuous(name="Wing Shape Integration (Ecc)")+
  scale_y_continuous(name="Proportion of Defects")+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population",
                      labels = c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))
s_2

# Correlation test for eccentricity vs proportion of defects
cor.test(mutant_integration_Eth$eccentricity, mutant_integration_Eth$defects_prop)
cor.test(mutant_integration_Zam$eccentricity, mutant_integration_Zam$defects_prop)


s_new <- ggplot(mutant_integration, aes( y = defects_prop, x = rSDE_alt , 
                                     color = population, shape = population)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.8, 0.9), 
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  scale_x_continuous(name="Wing Shape Integration (rSDE)")+
  scale_y_continuous(name="Proportion of Defects")+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population",
                      labels = c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))
s_new

s_2_new <- ggplot(mutant_integration, aes( y = defects_prop, x = rSDE2 , 
                                       color = population, shape = population)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.8, 0.9), 
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  scale_x_continuous(name="Wing Shape Integration (rSDE2)")+
  scale_y_continuous(name="Proportion of Defects")+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population",
                      labels = c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))
s_2_new

# Correlation test for rSDE vs proportion of defects
cor.test(mutant_integration_Eth$rSDE, mutant_integration_Eth$defects_prop)
cor.test(mutant_integration_Zam$rSDE, mutant_integration_Zam$defects_prop)

cor.test(mutant_integration_Eth$rSDE2, mutant_integration_Eth$defects_prop)
cor.test(mutant_integration_Zam$rSDE2, mutant_integration_Zam$defects_prop)


###############################################################

############# Using geomorph morphol.disparity function ########
wings_land<- wings[,10:105]
wings_descr <- wings[,c(3, 4, 5, 107)]
wings_csize <- wings$CS_rescaled

#converting landmarks into 3D array because geom
wings_land_3D <- arrayspecs(wings_land, 48, 2)

#converting to list, we have wings_land_3D from before
wings_list_3Darray <- list(wings_descr = wings_descr, 
                           wings_land_3D = wings_land_3D, 
                           wings_csize = wings_csize)

#converting to geomorph data frame
gdf <- geomorph.data.frame(shape = wings_land_3D,
                           population = wings_descr$population, 
                           sex = wings_descr$sex,
                           line = wings_descr$line,
                           size = wings_csize)



wings_allomAnova <- procD.lm(shape ~ size * sex * (population/line), 
                             data = gdf, iter = 100,
                             RRPP=TRUE, print.progress = FALSE)

######
wings_allomAnova2 <- procD.lm(shape ~ size * sex * population + (population/line), 
                             data = gdf, iter = 100,
                             RRPP=TRUE, print.progress = FALSE)


wings_allomAnova_nested2 <- anova(wings_allomAnova2, 
                                 error = c("Residuals", "population:line",
                                           "population:line", "population:line",
                                           "population:line", "population:line",
                                           "Residuals",
                                           "population:line")) 

#####
summary(wings_allomAnova)


wings_allomAnova_nested <- anova(wings_allomAnova, 
                                 error = c("Residuals", "Residuals",
                                  "population:line", "Residuals",
                                  "Residuals", "Residuals", "Residuals",
                                  "Residuals", "Residuals","Residuals",
                                  "Residuals")) 

MD <- morphol.disparity(wings_allomAnova2, groups = ~population, 
                        data = gdf, iter = 999)
summary(MD)



###############################################################

############## New Figure 3 Code ########################
legend_Figure3 <- get_legend(q_2 + theme(legend.position = "bottom"))

Figure3_new <- plot_grid(n, p , axis = "tblr", labels = "AUTO", label_size = 65)
Figure3_new_withLegend <- plot_grid(Figure3_new, legend_Figure3, ncol = 1, rel_heights = c(1, .1), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/Figure3_new_TotalVariance_rSDE.pdf", height = 19, width = 33)
Figure3_new_withLegend
dev.off()


########################################################

################# SupFigure 2 code ################

legend_SupFigure2 <- get_legend(s_2 + theme(legend.position = "bottom"))

SupFigure2 <- plot_grid(r, s , axis = "tblr", labels = "AUTO", label_size = 65)
SupFigure2_withLegend <- plot_grid(SupFigure2, legend_SupFigure2, ncol = 1, rel_heights = c(1, .1), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/SupFigure2_TotalVariance_Eccentricity_VS_Defects.pdf", height = 18, width = 32)
SupFigure2_withLegend
dev.off()

###################################################

################ SupFigure 5 ####################

legend_SupFigure2 <- get_legend(s_2 + theme(legend.position = "bottom"))

SupFigure5 <- plot_grid(s_new, s_2_new , axis = "tblr", labels = "AUTO", label_size = 65)
SupFigure5_withLegend <- plot_grid(SupFigure5, legend_SupFigure2, ncol = 1, rel_heights = c(1, .1), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/SupFigure5_rSDE_rSDE2_VS_Defects.pdf", height = 18, width = 32)
SupFigure5_withLegend
dev.off()

#################################################


################# SupFigure 6 code ################

SupFigure6 <- plot_grid(o, q , axis = "tblr", labels = "AUTO", label_size = 65)
SupFigure6_withLegend <- plot_grid(SupFigure6, legend_Figure3, ncol = 1, rel_heights = c(1, .1), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/SupFigure6_eccentricity_and_rSDE2.pdf", height = 19, width = 33)
SupFigure6_withLegend
dev.off()

###################################################
