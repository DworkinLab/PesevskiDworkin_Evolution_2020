######################################################################
###### Cell Density analysis and variation and Figure 4 script #######
######################################################################

########## Load Libraries and local functions #################
library(effects)
library(ggplot2)
library(lme4)
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

source('../misc/ID_LeveneStat_V1_2016.R', chdir = TRUE)
###############################################################


######################## Load Data ############################
#Size and shape data
source('./MP_African_ShortDataCleaning.R')

#Cell density data
source('./MP_African2016_CellCounts_DataCleanup.R')


###############################################################

#################### Set ggplot theme #########################

pref_theme <- theme_classic() + 
  theme(text = element_text(size = 65), 
        axis.title.y = element_text(margin = margin(t = 10, r = 20, b = 10, l = 10)),
        axis.title.x = element_text(margin = margin(t = 20, r = 10, b = 10, l = 10)),
        axis.text = element_text(margin = margin(t = 15, r = 15, b = 15, l = 15)))
theme_set(pref_theme)

###############################################################


################### Observation regions ######################

t <- image_read("../misc/Observation_positionsandregions_cropped.png")

u <- ggdraw() + draw_image(t, scale = 1)

###############################################################



################### Cell Density analysis across wing ######################

levels(cells$population) <- c("Ethiopia", "Zambia")

# Cell density mixed model with means for areas A, B, C, E, F and G. Trichome counts per 75x75 px box
lmm.1.AllAreas <- lmer(trichome_num ~ (wing_area + sex + population )^3 
                       + (1+sex|line:population) + (1|wing_ID) , 
                       data = cells)
summary(lmm.1.AllAreas)
Anova(lmm.1.AllAreas)

# Cell density mixed model using all measurement regions
lmm.2.AllAreas <- lmer(trichome_num ~ (wing_regions + sex + population )^3 
                       + (1+sex|line:population) + (1|wing_ID) , 
                       data = cells)
summary(lmm.2.AllAreas)


# Cell density mixed model using mean trichomes among all measurement reagions
lmm.1.Overall <- lmer(trichome_num ~ (sex + population )^3 
                       + (1+sex|line:population) + (1|wing_ID) , 
                       data = cells)
summary(lmm.1.Overall)
Anova(lmm.1.Overall)

lmm.1.Overall_marginal_means <- as.data.frame(Effect(c("sex", "population"), lmm.1.Overall))

#Supplementary Table 11
Anova(lmm.2.AllAreas)

# getting effects for plotting
wing_AllAreas1_marginal_means <- as.data.frame(Effect(c("wing_area", "sex", "population"), lmm.1.AllAreas))

wing_AllAreas2_marginal_means <- as.data.frame(Effect(c("wing_regions", "sex", "population"), lmm.2.AllAreas))


cells_lines <- aggregate(cells["trichome_num"], 
                         by = c(cells["sex"],
                                 cells["population"],
                                 cells["line"]),
                           mean, na.rm = FALSE)


## Figure 4B
v <- ggplot(wing_AllAreas1_marginal_means, aes(x=wing_area, y=fit, 
                                               shape = sex, colour=population))+
  geom_point(size = 14, position = position_dodge(0.9), show.legend = FALSE)+
  geom_errorbar(width=0.5, aes(ymin=lower, ymax=upper), 
                position = position_dodge(0.9), lwd = 1.5, show.legend = FALSE)+
  theme(legend.position = c(0.7, 0.15), 
        legend.box = "horizontal",
        legend.title = element_text(size = 38),
        legend.text = element_text(size = 36),
        legend.background = element_rect(fill="transparent")) +
  scale_y_continuous(name="Mean Cell Density")+
  scale_x_discrete(name="Wing Region")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population",
                      labels=c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name  ="Sex",
                       labels=c("Female", "Male"))

v

wing_AllAreas2_marginal_means$wing_regions<- droplevels(wing_AllAreas2_marginal_means$wing_regions)
levels(wing_AllAreas2_marginal_means$wing_regions) <- c("A1", "B1", "C1", "E1", "F1", "G1", "A2", "B2", "C2", "E2", "F2", "G2", "B3", "C3", "F3", "G3")


wing_AllAreas2_marginal_means$wing_regions <- factor(wing_AllAreas2_marginal_means$wing_regions, levels = sort(levels(wing_AllAreas2_marginal_means$wing_regions), decreasing = FALSE))


# Supplementary Figure S7
w <- ggplot(wing_AllAreas2_marginal_means, aes(x=wing_regions, y=fit, 
                                                 shape = sex, colour=population))+
  geom_point(size = 14, position = position_dodge(0.5))+
  geom_errorbar(width=0.5, aes(ymin=lower, ymax=upper), 
                position = position_dodge(0.5), lwd = 1)+
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        legend.text = element_text(size = 65),
        legend.background = element_rect(fill="transparent")) +
  scale_y_continuous(name="Mean Cell Density")+
  scale_x_discrete(name="Wing Region")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population",
                      labels=c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name  ="Sex",
                       labels=c("Female", "Male"))

w



v_alt <- ggplot(lmm.1.Overall_marginal_means, aes(x= population, y=fit, 
                                                   shape = sex, colour=population))+
  geom_point(size = 14, position = position_dodge(0.6), show.legend = FALSE)+
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), 
                position = position_dodge(0.6), lwd = 1.5, show.legend = FALSE)+
  geom_jitter(data = cells_lines, aes(x = population, y = trichome_num),
              size = 6, position = position_jitterdodge(jitter.width = 0.3), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        legend.title = element_text(size = 38),
        legend.text = element_text(size = 36),
        legend.background = element_rect(fill="transparent")) +
  scale_y_continuous(name="Cell Density")+
  scale_x_discrete(name="", labels = c("High-Altitude", "Low-Altitude"))+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population",
                      labels=c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name  ="Sex",
                       labels=c("Female", "Male"))

v_alt


###############################################################


###################### Cell Density CV ########################

### Cell density variability analysis

#Within individual variation, variation within each individual across the different regions

#because the mean cell density between regions is different I decided to do within individual CV to standardize for the differences in mean cell densities across the wing

#first I calculate the standard deviation and mean for each individual in each line and sex.

cells$individual<- as.factor(cells$individual)

cells_ind_sd<- as.data.frame(aggregate( cells[4], 
                                        c( cells["line"], 
                                           cells["sex"],
                                           cells["individual"]), 
                                        sd, na.rm=TRUE ))

names(cells_ind_sd)<- c("line", "sex", "individual", "sd_trichome")

cells_ind_mean<- as.data.frame(aggregate(cells[4],
                                         c( cells["line"], 
                                            cells["sex"],
                                            cells["individual"]), 
                                         mean, na.rm=TRUE ))

names(cells_ind_mean)<- c("line", "sex", "individual", "mean_trichome")

# Calculating the within individual CV by dividing the individual standard deviation by the individual mean.
cells_combined <- merge(cells_ind_sd,cells_ind_mean)

cells_combined$cellsCV <- cells_combined$sd_trichome/cells_combined$mean_trichome

#adding a population factor
pop_dummy <- substr(cells_combined$line, 1,1)
cells_combined$population <- ifelse(pop_dummy == "Z", "Zambia", "Ethiopia")
cells_combined$population <- factor(cells_combined$population)

#calculating the line means for within individual CV
cells_lineCVmeans<- as.data.frame(aggregate( cells_combined[6],  
                                             c( cells_combined["line"], 
                                                cells_combined["sex"]), 
                                             mean, na.rm=TRUE ))

#adding a population factor
pop_dummy <- substr(cells_lineCVmeans$line, 1,1)
cells_lineCVmeans$population <- ifelse(pop_dummy == "Z", "Zambia", "Ethiopia")
cells_lineCVmeans$population <- factor(cells_lineCVmeans$population)


# Cell density CV model for within individual CV
lmm.CVCells_gamma <- glmer(cellsCV ~ population*sex + (1 + sex|line),
                          data = cells_combined, family = "Gamma",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))

summary(lmm.CVCells_gamma)

# Supplementary Table S12
Anova(lmm.CVCells_gamma)


lmm.CVCells_gamma_marginal_means <- as.data.frame(Effect(c("population", "sex"), 
                                                         lmm.CVCells_gamma))

x <- ggplot(lmm.CVCells_gamma_marginal_means, aes( y = fit, x = population, 
                                    colour=population, shape = sex)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = cells_lineCVmeans, aes(x = population, y = cellsCV),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Cell Density Variation (CV)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Sex", 
                       labels = c("Female", "Male"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
x



# Now to check if its consistent we will first average cell density by line and then calculate within line CV

cells_line_sd<- as.data.frame(aggregate( cells[4], 
                                        c( cells["line"], 
                                           cells["sex"]), 
                                        sd, na.rm=TRUE ))

names(cells_line_sd)<- c("line", "sex", "sd_trichome")

cells_line_mean<- as.data.frame(aggregate(cells[4],
                                         c( cells["line"], 
                                            cells["sex"]), 
                                         mean, na.rm=TRUE ))

names(cells_line_mean)<- c("line", "sex", "mean_trichome")

#Calculating the within line CV by dividing the line standard deviation by the line mean.
cells_lines_combined <- merge(cells_line_sd, cells_line_mean, by = c("line", "sex"))

cells_lines_combined$cellsLinesCV <- cells_lines_combined$sd_trichome / cells_lines_combined$mean_trichome

#adding a population factor
pop_dummy <- substr(cells_lines_combined$line, 1,1)
cells_lines_combined$population <- ifelse(pop_dummy == "Z", "Zambia", "Ethiopia")
cells_lines_combined$population <- factor(cells_lines_combined$population)


# checking how correlated the two ways to calculate cv are
CVind_CVline <- merge(cells_lines_combined, cells_lineCVmeans, by = c("line", "sex", "population"))

check <- ggplot(data=CVind_CVline, aes(y = cellsCV, x = cellsLinesCV, 
                                       color = population)) +
  geom_point(size = 10, show.legend = TRUE)  +
  theme(legend.position = c(0.3, 0.9))+
  geom_smooth(method='lm',formula=y~x)+
  scale_x_continuous(name="Cell densuty CV") +
  scale_y_continuous(name = "Cell Density Lev Dev")+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population",
                      labels = c("High-Altitude", "Low-Altitude"))

check

cor.test(CVind_CVline$cellsLinesCV, CVind_CVline$cellsCV) #r = 0.91

lmm.CV.linesCells_gamma <- glmer(cellsLinesCV ~ population*sex + (1|line:population),
                           data = cells_lines_combined, family = "Gamma",
                           control=glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun=2e5)))

summary(lmm.CV.linesCells_gamma)
Anova(lmm.CV.linesCells_gamma)
plot(allEffects(lmm.CV.linesCells_gamma), multiline = TRUE)

lmm.CVlinesCells_gamma_marginal_means <- as.data.frame(Effect(c("population", "sex"), 
                                                              lmm.CV.linesCells_gamma))


y <- ggplot(lmm.CVlinesCells_gamma_marginal_means, aes( y = fit, x = population, 
                                                   colour=population, shape = sex)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = cells_lines_combined, aes(x = population, y = cellsLinesCV),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Cell Density Variation (CV)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Sex", 
                       labels = c("Female", "Male"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
y


#small differences between CV calculated by ind first and then averaged by line and cv calculated by line. the overall pattern is consistent, There is a sex effect but no population effect.

###############################################################


####################### Levene Deviates #######################
# Calculating Within line variation for cell density using Levene Deviats
cells$LevDev <- with(cells, LeveneDeviates(trichome_num, 
                                           group = interaction(line, sex, wing_ID, drop=TRUE),
                                           med  = TRUE, log_trans = TRUE))


cells$LevDev2 <- cells$LevDev + 0.0001

lmm.cell_LevDev_gamma <- glmer(LevDev2 ~ sex*population + (1|line:population),
                          data = cells, family = "Gamma",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))

summary(lmm.cell_LevDev_gamma)

# Supplementary Table S13
Anova(lmm.cell_LevDev_gamma)

cells_LineLevDev <- aggregate( cells[,c(4,5,11, 12, 13, 14, 15)], 
                               c(cells["sex"],
                                 cells["line"], 
                                 cells["population"]), 
                               mean, na.rm=TRUE )

lmm.cell_LevDev_gamma_MargMean <- as.data.frame(Effect(c("sex", "population"), lmm.cell_LevDev_gamma ))


# Supplementary Figure S8A
z <- ggplot(lmm.cell_LevDev_gamma_MargMean, aes( y = fit, x = population, 
                                             colour=population, shape = sex)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=lower, ymax=upper), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = cells_LineLevDev, aes(x = population, y = LevDev2),
              size = 10, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y=" Cell Density variation (LD)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Sex", 
                       labels = c("Female", "Male"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
z

###############################################################

################# Cell Density LevDev vs CV ###################
# Checking how correlated CV and Lev Dev are
Cells_CV_LevDev <- merge(cells_lines_combined, cells_LineLevDev)

aa <- ggplot(data=Cells_CV_LevDev, aes(y = LevDev, x = cellsLinesCV, 
                                      color = population)) +
  geom_point(size = 10, show.legend = FALSE)  +
  theme(legend.position = c(0.8, 0.9))+
  scale_x_continuous(name="Cell Density Variation (CV)") +
  scale_y_continuous(name = "Cell Density Variation (LD)")+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", "Low-Altitude"))

aa

###############################################################

########### Cell density CV vs Wing size CV, Wing Shape Total Variance and Eccentricty, mutant frequency ###################

## Cell Density CV vs Wing size CV
source('./MP_AfricanCanalization_Figure2_WingSize_CV_LevDev_WingDefects.R')

cv_lines[,c(1, 2, 6, 7)] <- lapply(cv_lines[,c(1, 2, 6, 7)], toupper)

CVsize_CVcells <- merge(cells_lineCVmeans, cv_lines, by = c("sex", "line", "population"))


CVsize_CVcells$cv <- as.numeric(CVsize_CVcells$cv)
#plotting line CVs for wing size against line CVs for cell density
#plot(CVsize_CVcells$cv,CVsize_CVcells$cellsCV)

###Figure 4D####
ab <- ggplot(CVsize_CVcells, aes( y = cellsCV, x = cv, 
                                 color = population)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = c(0.6, 0.9), legend.background = element_rect(fill="transparent"))+
  labs(x="Wing Size Variation (CV)", y="Cell Density Variation (CV)")+
  scale_color_manual(name= "Population",
                     values = c("black","gray60"), 
                     labels = c("High-Altitude", "Low-Altitude")) 

ab

CVsize_CVcells_eth <- subset(CVsize_CVcells, population == "Ethiopia")

CVsize_CVcells_zam <- subset(CVsize_CVcells, population == "Zambia")

cor.test(CVsize_CVcells_eth$cv, CVsize_CVcells_eth$cellsCV)
cor.test(CVsize_CVcells_zam$cv, CVsize_CVcells_zam$cellsCV)

## Cell Density CV vs Wing Shape Total variance and Eccentricty

source('./MP_AfricanCanalization_Figure3_WingShape_TotVar_Ecc_Defects.R')

integ_CVcells <- merge(cells_lineCVmeans,integration_measures2)

integ_CVcells_Eth <- subset(integ_CVcells, population == "Ethiopia")

integ_CVcells_Zam <- subset(integ_CVcells, population == "Zambia")


# Supplementary Figure S8C
ac <- ggplot(integ_CVcells, aes( y = cellsCV, x = tot_var_alt, 
                                color = population)) +
  geom_point(size=14, show.legend = FALSE)+
  theme(legend.position = c(0.8, 0.9), 
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  labs(x="Wing Shape Variation (TV)", y="Cell Density Variation (CV)")+
  scale_color_manual(name= "Population", 
                     labels =c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))+
  scale_x_continuous(limits = c(0.11, 0.26))
ac


cor.test(integ_CVcells_Eth$cellsCV, integ_CVcells_Eth$tot_var_alt)
cor.test(integ_CVcells_Zam$cellsCV, integ_CVcells_Zam$tot_var_alt)


# Supplementary Figure S8D
ad <- ggplot(integ_CVcells, aes( y = cellsCV, x = eccentricity, color = population)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = c(0.8, 0.9), 
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  labs(x="Wing Shape Integration (Ecc)", y="Cell Density Variation (CV)")+
  scale_color_manual(name= "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))

ad

cor.test(integ_CVcells_Eth$cellsCV, integ_CVcells_Eth$eccentricity)
cor.test(integ_CVcells_Zam$cellsCV, integ_CVcells_Zam$eccentricity)

# #Cell density CV vs rSDE
# 
# ad_2 <- ggplot(integ_CVcells, aes( y = cellsCV, x = rSDE_alt, color = population)) +
#   geom_point(size = 14, show.legend = FALSE)+
#   theme(legend.position = c(0.8, 0.9), 
#         legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
#   labs(x="rSDE", y="Cell density CV")+
#   scale_color_manual(name= "Population", 
#                      labels = c("High-Altitude","Low-Altitude"), 
#                      values=c("black","gray60"))+
#   scale_shape_discrete(name  ="Population", 
#                        labels = c("High-Altitude", "Low-Altitude"))
# 
# ad_2
# 
# cor.test(integ_CVcells_Eth$cellsCV, integ_CVcells_Eth$rSDE_alt)
# cor.test(integ_CVcells_Zam$cellsCV, integ_CVcells_Zam$rSDE_alt)

#############################################################

###### Cell Density CV vs Proportion of defects#########
mutant_frequencies_averaged <- aggregate(x = mutant_frequencies[,c(11, 12, 13 )], 
                                         by = c(mutant_frequencies["population"], mutant_frequencies["line"]), 
                                         FUN = mean )


mutant_cellCV <- merge(x = cells_lineCVmeans, y = mutant_frequencies_averaged )


mutant_cellCV_Eth <- subset(mutant_cellCV, population == "Ethiopia")

mutant_cellCV_Zam <- subset(mutant_cellCV, population == "Zambia")

# Supplementary Figure S8B
ae <- ggplot(mutant_cellCV, aes( y = cellsCV, x = defects_prop, 
                                 color = population)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = c(0.8, 0.9), 
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  labs(y= "Cell Density Variation (CV)", x="Proportion of Defects")+
  scale_color_manual(name= "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  ="Population", 
                       labels = c("High-Altitude", "Low-Altitude"))
ae


cor.test(mutant_cellCV_Eth$cellsCV, mutant_cellCV_Eth$defects_prop)
cor.test(mutant_cellCV_Zam$cellsCV, mutant_cellCV_Zam$defects_prop)

################################################################################

################# Figure 4 code ############
legend_Figure4 <- get_legend(w + theme(legend.position = "bottom", legend.box = "vertical"))

Figure4 <- plot_grid(u, v_alt, x, ab, axis = "tblr", labels = "AUTO", label_size = 65 )

Figure4_withLegend <- plot_grid(Figure4, legend_Figure4, ncol = 1, rel_heights = c(1, .1), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/Figure4_CellDensity_Per_Region_and_CellDensity_CV.pdf", height = 35,  width = 38)
Figure4_withLegend
dev.off()
############################################


############# Supplemental Figure 7 ########

SupFigure7 <- plot_grid(w, nrow = 1, axis = "tblr")


pdf(file = "../outputs/PaperFigures/SupFigure7_CellDensity_allMeasurementRegions_DifferencesInCellDensity.pdf", height = 16,  width = 32)
SupFigure7
dev.off()

###########################################

########## Supplemental Figure 8 ##########


legend_SupFigure8 <- get_legend(w + theme(legend.position = "bottom", legend.box = "vertical"))

SupFigure8 <- plot_grid(z, ae, ac, ad, axis = "tblr", labels = "AUTO", label_size = 65)
SupFigure8_withLegend <- plot_grid(SupFigure8, legend_SupFigure8, ncol = 1, rel_heights = c(1, .1), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/SupFigure8_TotalVariance_Eccentricity_VS_CellDensityCV.pdf", height = 35, width = 38)
SupFigure8_withLegend
dev.off()

###########################################

