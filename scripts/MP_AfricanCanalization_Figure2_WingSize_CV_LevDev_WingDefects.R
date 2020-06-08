########################################################################
### Wing size Integration Measures, Wing defects and Figure 2, Supplementary Figures 1, 3 and 4 script ####
########################################################################

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

#### Microenvironmental Canalization - Wing size CV #####

##first we center the wing size
df3.mod.size <- lm(CS_rescaled ~ sex*population, data = wings)
wings$centered.size <- coef(df3.mod.size)[1] + df3.mod.size$resid

# CV function
cv <- function(x = CS_rescaled) {sd(x)/mean(x)}

# Calculate CV for each line
cv_lines <- with(wings, 
                 tapply(CS_rescaled, INDEX=interaction(population,line, sex, drop=TRUE), cv))

# Average wing size for each line
mean_lines <- with(wings, 
                   tapply(CS_rescaled, INDEX=interaction(population,line, sex, drop=TRUE), mean))

# standard deviation of each line
sd_lines <- with(wings, 
                 tapply(CS_rescaled, INDEX=interaction(population,line, sex, drop=TRUE), sd))

# Wing size CV data frame and clean up 
cv_lines <- data.frame(as.character(row.names(cv_lines)), cv_lines, mean_lines, sd_lines)
colnames(cv_lines)[1] <- "pop.line.sex"

cv_lines$population <- matrix(unlist(
  strsplit(as.character(row.names(cv_lines)), split = '.', fixed =TRUE)), 
  ncol =3, byrow=TRUE)[,1]
cv_lines$population <- factor(cv_lines$population)


cv_lines$sex <- matrix(unlist(
  strsplit(as.character(row.names(cv_lines)), split = '.', fixed =TRUE)), 
  ncol =3, byrow=TRUE)[,3]
cv_lines$sex <- factor(cv_lines$sex)

cv_lines$line <- matrix(unlist(
  strsplit(as.character(row.names(cv_lines)), split = '.', fixed =TRUE)), 
  ncol =3, byrow=TRUE)[,2]
cv_lines$line <- factor(cv_lines$line)

# Are the Ethiopian lines more variable?
cv_populations <- aggregate(cv_lines[,2:4],
                            by = c(cv_lines["population"],
                                   cv_lines["sex"]),
                            mean, na.rm = FALSE)

# check the relationship between mean and sd
with(cv_lines[cv_lines$sex == "m",], plot(mean_lines, sd_lines))


### CV glmer model
lmm.CV_gamma <- glmer(cv_lines ~ population*sex + (1|line:population),
                              data = cv_lines, family = "Gamma",
                              control=glmerControl(optimizer="bobyqa",
                                                   optCtrl=list(maxfun=2e5)))
summary(lmm.CV_gamma)


# Table S6
Anova(lmm.CV_gamma)

#Extracting effects
lmm.CV_gamma_marginal_means <- as.data.frame(Effect(c("population", "sex"), 
                                                    lmm.CV_gamma))


# Figure 2A: CV by sex and pop
f <- ggplot(lmm.CV_gamma_marginal_means, aes( y = fit, x = population, 
                                                        colour=population, shape = sex)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=lower, ymax=upper), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = cv_lines, aes(x = population, y = cv_lines),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Wing Size Variation (CV)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Sex", 
                       labels = c("Female", "Male"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
f


# for legend
f_1 <- ggplot(lmm.CV_gamma_marginal_means, aes( y = fit, x = population, 
                                                colour=population, shape = sex)) +
  geom_point(size = 14, show.legend = TRUE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=lower, ymax=upper), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = cv_lines, aes(x = population, y = cv_lines),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Wing Size Variation (CV)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Sex", 
                       labels = c("Female", "Male"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
f_1


###############################################################

################### Wing size LevDev ##########################

# Generate the Levene's Deviates by line and sex. We included sex, because even after mean centering for the sexual dimorphism, there was still evidence of additional variance in females. Likely this is due to line specific SSD that has to be modeled out with the random effects. 
wings$LevDev <- with(wings, LeveneDeviates(CS_rescaled, 
                                           group = interaction(line, sex, drop=TRUE),
                                           med  = TRUE, log_trans = TRUE))

wings$LevDev2 <- wings$LevDev + 0.0001


# Levene Deviates glmer model
lmm.LevDev_gamma <- glmer(LevDev2 ~ sex*population + (1|line:population),
                      data = wings, family = "Gamma",
                      control=glmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=2e5)))

summary(lmm.LevDev_gamma)

# Table S7
Anova(lmm.LevDev_gamma)

# extracting effects
lmm.LevDev_gamma_MargMean <- as.data.frame(Effect(c("sex", "population"), lmm.LevDev_gamma ))

# averaging lev dev by line
wings_LineLevDev <- aggregate( wings[,109:112], 
                               c(wings["sex"],
                                 wings["line"], 
                                 wings["population"]), 
                               mean, na.rm=TRUE )


# Supplementary Figure 3A: Levene Deviates
g <- ggplot(lmm.LevDev_gamma_MargMean, aes( y = fit, x = population, 
                                                    colour=population, shape = sex)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=lower, ymax=upper), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = wings_LineLevDev, aes(x = population, y = LevDev2),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Wing size Variation (LD)")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Sex", 
                       labels = c("Female", "Male"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
g

###############################################################

######################## Wing defects #########################

mutant_frequencies <- read.csv("../data/African_Mutation_frequencies.csv", h = T)

mutant_frequencies$alv_prop <- with(mutant_frequencies, alv_num/total_wings)
mutant_frequencies$micv_prop <- with(mutant_frequencies, micv_num/total_wings)
mutant_frequencies$defects_prop <- with(mutant_frequencies, total_defects/total_wings)

# create a population variable 
boolean_for_line <- grepl(pattern = "EF", mutant_frequencies$line)
mutant_frequencies$population <- as.factor(ifelse(boolean_for_line, "Ethiopia", "Zambia"))

# Wing Defects glm
mut.glm1 <- glm(cbind(total_defects, total_wings) ~ population*sex, data = mutant_frequencies , family = binomial) 
summary (mut.glm1)

# Table S8
Anova(mut.glm1)

# Wing defects effects
Defects_marg_means <- as.data.frame(Effect(c("population", "sex"), mut.glm1))

# Supplementary Figure S1B
h <- ggplot(Defects_marg_means, aes( y = fit, x = population, 
                                              colour=population, shape = sex)) +
  geom_point(size = 14, show.legend = FALSE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=lower, ymax=upper), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = mutant_frequencies, aes(x = population, y = defects_prop),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "horizontal", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Proportion of Defects")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Sex", 
                       labels = c("Female", "Male"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
h


h_2 <- ggplot(Defects_marg_means, aes( y = fit, x = population, 
                                       colour=population, shape = sex)) +
  geom_point(size = 14, show.legend = TRUE,  position = position_dodge(0.8)) +
  geom_errorbar(width=0.2, aes(ymin=lower, ymax=upper), lwd = 1.5, 
                position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = mutant_frequencies, aes(x = population, y = defects_prop),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 36),
        legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  labs(x="", y="Proportion of Defects")+
  scale_color_manual(name = "Population",
                     labels = c("High-Altitude", "Low-Altitude"),
                     values=c("black","gray60")) +
  scale_shape_discrete(name  ="Sex", 
                       labels = c("Female", "Male"))+
  scale_x_discrete(labels=c("High-Altitude", "Low-Altitude"))
h_2

# Supplementary Figure S1A
# loading cropped version of wingplot as a pdf
i_1 <- image_read("../outputs/Example_WingDefects.png")

# turing wingplot pdf into ggplot friendly object
i <- ggdraw() + draw_image(i_1)

###############################################################


################# Wing size cv vs defects #####################
# Clean up to merge CV data frame and defects data frame
cv_lines$sex <- toupper(cv_lines$sex)

cv_lines$line <- matrix(unlist(
  strsplit(as.character(row.names(cv_lines)), split = '.', fixed =TRUE)), 
  ncol =3, byrow=TRUE)[,2]

cv_lines$line <- as.factor(toupper(cv_lines$line))
names(cv_lines)[2] <- "cv"

cv_lines_2 <- cv_lines[, c("line", "sex", "cv")]

# merge cvs and mutant frequencies
mutant_CV <- data.frame(merge(x = cv_lines_2, y = mutant_frequencies ) )


#Figure2B
j <- ggplot(data=mutant_CV, aes(y = cv, x = defects_prop, 
                               color = population)) +
  geom_point(size = 10, show.legend = FALSE)  +
  theme(legend.position = c(0.8, 0.9))+
  scale_x_continuous(name="Proportion of Defects") +
  scale_y_continuous(name = "Wing Size Variation (CV)", 
                     limits = c(0.01, 0.041))+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", "Low-Altitude"))
 
j

# Correlation Tests for wing size CV vs proportion of defects for each population
mutant_CV_eth <- subset(mutant_CV, population == "Ethiopia")

cor.test(mutant_CV_eth$cv, mutant_CV_eth$defects_prop)

mutant_CV_zam <- subset(mutant_CV, population == "Zambia")

cor.test(mutant_CV_zam$cv, mutant_CV_zam$defects_prop)

###############################################################

################# Wing size Lev Dev vs defects #####################

# Correlation beteeen Proportion of defects and wing defects
wings_LineLevDev$sex <- toupper(wings_LineLevDev$sex)
wings_LineLevDev$line <- toupper(wings_LineLevDev$line)

mutant_levDev <- data.frame(merge(wings_LineLevDev, mutant_CV,  by = c("sex", "population", "line")))

# Supplementary Figure S4
k <- ggplot(data=mutant_levDev, aes(y = LevDev, x = defects_prop, 
                               color = population, shape = population)) +
  geom_point(size = 10, show.legend = TRUE)  +
  theme(legend.position = "bottom")+
  scale_x_continuous(name="Proportion of Defects") +
  scale_y_continuous(name = "Wing Size Variation (LD)")+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", "Low-Altitude")) +
  scale_shape_discrete(name  ="Population", 
                      labels = c("High-Altitude", "Low-Altitude"))
  

k

# Correlation test for wing size Levene Deviates vs poportion of defects
mutant_levDev_eth <- subset(mutant_levDev, population == "Ethiopia")

cor.test(mutant_levDev_eth$LevDev, mutant_levDev_eth$defects_prop)

mutant_levDev_zam <- subset(mutant_levDev, population == "Zambia")

cor.test(mutant_levDev_zam$LevDev, mutant_levDev_zam$defects_prop)
#################################################################


########## Correlation between CV and Lev Dev ###################
# check to see how correlated CV and Lev Dev are
l <- ggplot(data=mutant_levDev, aes(y = LevDev, x = cv, 
                                      color = population)) +
  geom_point(size = 10, show.legend = FALSE)  +
  theme(legend.position = c(0.8, 0.9))+
  scale_x_continuous(name="Wing Size Variation (CV)") +
  scale_y_continuous(name = "Wing Size Variation (LD)")+
  scale_colour_manual(values = c("black","gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", "Low-Altitude"))

l

# correlation test between wing size CV and Lev Dev
cor.test(mutant_levDev_eth$LevDev, mutant_levDev_eth$cv)
cor.test(mutant_levDev_zam$LevDev, mutant_levDev_zam$cv)

###############################################################

####################### Wing Size GCV -  no longer included #########################

# Calculating the Genetic Coefficient of Variation for each population

#Removing Sex variation
lmm.size.zambia <- lmer(CS_rescaled ~ sex + (1 + sex|line), 
                        data = wings, subset = population == "Zambia")


lmm.size.ethiopia <- lmer(CS_rescaled ~ sex + (1 + sex|line), 
                          data = wings, subset = population == "Ethiopia")

# Coefficient of genetic variation for size for Zambia is:
CGV.zam<-as.data.frame(VarCorr(lmm.size.zambia))[1,5]/fixef(lmm.size.zambia)[1]
# Coefficient of genetic variation for size for Ethiopia is:
CGV.eth<-as.data.frame(VarCorr(lmm.size.ethiopia))[1,5]/fixef(lmm.size.ethiopia)[1]

# Making a dataframe for CGV values
CGV<-as.data.frame(cbind(Zambia = CGV.zam, Ethiopia = CGV.eth))
CGV_row <- as.data.frame(rbind(Zambia = CGV.zam, Ethiopia = CGV.eth))

#Calculating CIs for CGV for each population
CIs_zambia <- confint(lmm.size.zambia)
CIs_ethiopia <- confint(lmm.size.ethiopia)

# # Approximate lower bounds for zambia
zam.low.CI<-CIs_zambia[1,1]/CIs_zambia[5,1]
# Approximate upper bounds.
zam.up.CI<- CIs_zambia[1,2]/CIs_zambia[5,2]

zam.CI<-as.data.frame(rbind(LowerBound=zam.low.CI,UpperBound=zam.up.CI))
colnames(zam.CI)<-"Zambia"

zam.ci.rows <- as.data.frame(cbind(LowerBound=zam.low.CI,UpperBound=zam.up.CI))

# # Approximate lower bounds for Ethiopia
eth.low.CI<-CIs_ethiopia[1,1]/CIs_ethiopia[5,1]
# 
# # Approximate upper bounds.
eth.up.CI<-CIs_ethiopia[1,2]/CIs_ethiopia[5,2]

eth.CI<-as.data.frame(rbind(LowerBound=eth.low.CI,UpperBound=eth.up.CI))
colnames(eth.CI)<-"Ethiopia"

eth.CI.rows<-as.data.frame(cbind(LowerBound=eth.low.CI,UpperBound=eth.up.CI))


Afr.CIs<-cbind(zam.CI,eth.CI)
Afr.CIs.rows <- rbind(zam.ci.rows,eth.CI.rows)

CGV_CIs <- as.data.frame(cbind(CGV_row, Afr.CIs.rows))

CGV_CIs$population <- c("Zambia", "Ethiopia")


names(CGV_CIs) <- c("cgv","LowerBound", "UpperBound", "population")


# CGV by population - not included in the paper
m <- ggplot(data = CGV_CIs, aes(y = cgv, x = population, 
                                colour = population)) +
  geom_point(size = 14, show.legend = FALSE) +
  geom_errorbar(aes(ymin=LowerBound, ymax=UpperBound), width = 0.2, lwd = 1.2, show.legend = FALSE) +
  theme(legend.position = c(0.5, 0.9)) +
  labs(y = "CGV", x = "") +
  scale_x_discrete(name="", 
                   labels = c("High-Altitude","Low-Altitude")) +
  scale_colour_manual(name = "Population",
                      labels = c("High-Altitude","Low-Altitude"),
                      values = c("black","gray60")) +
  scale_shape_discrete(name = "Population",
                       labels = c("High-Altitude","Low-Altitude"))
m
###############################################################


################# Figure 2 ##################
#Figure 2
legend_Figure2 <- get_legend(f_1 + theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 65), legend.title = element_text(size = 65)))

Figure2 <- plot_grid(f,j, align = "hv", labels = "AUTO", ncol = 2,  label_size = 65)

Figure2_withLegend <- plot_grid(Figure2, legend_Figure2, ncol = 1, rel_heights = c(1, .2), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/Figure2_Coefficient_of_variation_wingsize_and_proportion_of_defects.pdf", height = 18, width = 32)
Figure2_withLegend
dev.off()
###################################################


################# Supplementary Figure 1 ##################

SupFig1 <- get_legend(h_2 + theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 65), legend.title = element_text(size = 65)))

SupFigure1 <- plot_grid(i, h_2, labels = "AUTO",  label_size = 65, rel_widths = c(1, 1))

pdf(file = "../outputs/PaperFigures/SupFigure1_WingDefect_Examples_and_By_Population.pdf", height = 15, width = 30)
SupFigure1
dev.off()

################################################################

################# Supplementary Figure 3 ##################

SupFigure3 <- plot_grid(g, l, labels = "AUTO",  label_size = 65, rel_widths = c(1, 1))

SupFigure3_withLegend <- plot_grid(SupFigure3, ncol = 1, legend_Figure2, rel_heights = c(1, 0.2), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/SupFigure3_WingSize_LeveneStatistic_by_pop_LevDev_vs_CV.pdf", height = 18, width = 32)
SupFigure3_withLegend
dev.off()
################################################################

################# Supplementary Figure 4 ##################

SupFigure4 <- k

pdf(file = "../outputs/PaperFigures/SupFigure4_LevDev_vs_Defects_and_CoefGenVar.pdf", height = 15, width = 20)
SupFigure4
dev.off()
################################################################