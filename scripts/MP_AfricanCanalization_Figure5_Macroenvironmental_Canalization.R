######################################################################
########## MacroEnvironmental Canalization Figure 5 script ###########
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

source('../misc/WINGPLOTSOURCE.R', chdir = TRUE)

###############################################################

######################## Load Data ############################
#Size and shape data
source('./MP_MacroEnvironmental_Canalization_Cleanup.R')

###############################################################

#################### Set ggplot theme #########################

pref_theme <- theme_classic() + 
  theme(text = element_text(size = 65), 
        axis.title.y = element_text(margin = margin(t = 10, r = 20, b = 10, l = 10)),
        axis.title.x = element_text(margin = margin(t = 20, r = 10, b = 10, l = 10)),
        axis.text = element_text(margin = margin(t = 15, r = 15, b = 15, l = 15)))
theme_set(pref_theme)

###############################################################

##################### Pop means ######################

wings_popMeans <- aggregate( wings_temp[,c(7:102,103)], 
                             c(wings_temp["sex"],
                               wings_temp["population"], 
                               wings_temp["temperature"],
                               wings_temp["temp_num"],
                               wings_temp["temp_zeroed"]), 
                             mean, na.rm=TRUE )

###############################################################

##################### Removing Lines with low Ns ###############
with(wings_temp, table(line, temperature))

#Getting rid of lines:temperature that have less than 20 individuals (males and females together)

LowN_index <- with(wings_temp, line:temperature %in% c("e117:t18", "e117:t28", 
                                                       "e122:t18",  "e19:t18", 
                                                       "e19:t28", "e65:t28", 
                                                       "z186:t24", "z216:t28",  
                                                       "z360:t28",  "z366:t18", 
                                                       "z366:t24", "z366:t28",  
                                                       "z367:t18"))

wings_temp_NoLowNs <- wings_temp[(!LowN_index), ]

with(wings_temp_NoLowNs, table(line, temperature:sex))



wings_temp <- wings_temp_NoLowNs
wings_temp$line <- droplevels(wings_temp$line)

with(wings_temp, table(line, temperature:sex))

################################################################


################## Wing size reaction norms ###################
wings_temp$population <- relevel(wings_temp$population, "Zambia")

temp_mod1 <- lmer(CS ~ poly(temp_zeroed, 2)*sex*population 
                  + (1 + poly(temp_zeroed,2) + sex|population:line), 
                  data = wings_temp)

summary(temp_mod1)

# Supplementary Table S14
car::Anova(temp_mod1)

ReactionNorms_marginal_means <- as.data.frame(Effect(c("temp_zeroed", "sex", "population"), temp_mod1))

# Figure 5A
af <- ggplot(data = ReactionNorms_marginal_means, 
             aes(x = (temp_zeroed + 18), y = fit, colour = population, shape = population, linetype = sex)) +
  geom_line(aes(group = population:sex), lwd = 2, show.legend = TRUE) +
  geom_ribbon(aes(ymin=lower,ymax=upper, group = population:sex), alpha=0.15, 
              linetype = 0, show.legend = FALSE) +
  geom_point(data = wings_popMeans, aes(x = (temp_zeroed + 18), 
                                        y = CS, colour = population), size = 14) +
  theme(legend.position = c(0.3, 0.1), legend.box = "horizontal", 
        legend.key.width = unit(1.5, "cm"), 
        legend.background = element_rect(fill = alpha("white", 0), size=0.5),
        legend.title = element_text(size = 38),
        legend.text = element_text(size = 36))+
  scale_x_continuous(name = expression("Temperature " ( degree*C)) , breaks = c(18, 24, 28)) +
  ylab("Wing Size")+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name = "Population",
                       labels = c("High-Altitude", "Low-Altitude"))+
  scale_linetype_discrete(name = "Sex",
                          labels = c("Female", "Male"))
af

###############################################################


################ Wing shape Reaction norms#####################
wings_land<- wings_temp[,7:102]
wings_descr <- wings_temp[,c(105:114)]
wings_csize <- wings_temp$CS

#converting landmarks into 3D array because geom
wings_land_3D <- arrayspecs(wings_land, 48, 2)


gdf <- geomorph.data.frame(shape = wings_land_3D,
                           population = wings_descr$population,
                           sex = wings_descr$sex,
                           line = wings_descr$line,
                           size = wings_csize,
                           nutrition = wings_descr$nutrition,
                           temperature = wings_descr$temperature,
                           temp_zeroed = wings_descr$temp_zeroed,
                           temp_zeroed_sq = (wings_descr$temp_zeroed)^2)




#full model
wings_allomAnova <- procD.lm(shape ~ (size + sex + temp_zeroed + population)^3 + temp_zeroed_sq 
                             + size:temp_zeroed_sq + sex:temp_zeroed_sq + population:temp_zeroed_sq 
                             + size:sex:temp_zeroed_sq + size:population:temp_zeroed_sq 
                             + sex:population:temp_zeroed_sq +
                             + population/line,
                             data = gdf, iter = 2000,
                             RRPP=TRUE, print.progress = TRUE)

summary(wings_allomAnova)

# updating the residuals to account for nested design
wings_allomAnova_nested <- anova(wings_allomAnova, 
                                 error = c("Residuals", "Residuals", "Residuals",
                                           "population:line", "Residuals",
                                           "Residuals", "Residuals",
                                           "Residuals","Residuals", 
                                           "Residuals","Residuals",
                                           "Residuals", "Residuals",
                                           "Residuals", "Residuals",
                                           "Residuals", "Residuals",
                                           "Residuals", "Residuals",
                                           "Residuals", "Residuals",
                                           "Residuals")) 

# Supplementary Table S17
summary(wings_allomAnova_nested)

# NOTE: The output from this model is slightly different from the output reported in the manuscript because that MANOVA table was produced before we removed lines at different temperatures that had low Ns. 
###############################################################


########################## Temp PCA ###########################

wing_temp_PCs <- prcomp(wings_temp[,7:102])

wingstempPC <- data.frame(wings_temp, wing_temp_PCs$x[,1:58])

###############################################################.


################ Temperature Wing Mutations ###################
NewMutant <- read.csv("../data/African_Condition_2_BinaryFormatMutants_macrocanalization.csv")

# cleaning frequency of wing defects data for wings raised at different temperatures
# fixing predictor variables
NewMutant$line <- as.factor(NewMutant$line)
NewMutant$condition <- as.factor(NewMutant$condition)
NewMutant$sex <- as.factor(NewMutant$sex)

levels(NewMutant$condition) <- c("f100", "f15", "f5", "f5","t18", "t24", "t28")

# creating a nutrution variable
nutrition <- NewMutant$condition

levels(nutrition) <- c("f100", "f15", "f5","f100", "f100", "f100")

# creatring a temperarture variable
temperature <- NewMutant$condition

levels(temperature) <- c("t24", "t24", "t24","t18", "t24", "t28")

# adding nutrution and temperature variables
NewMutant <- cbind(NewMutant, nutrition, temperature)

NewMutant$temperature <- relevel(NewMutant$temperature, "t18")

# adding a population variable
pop_dummy <- substr(NewMutant$line, 1, 1)

NewMutant$population <- ifelse(pop_dummy == "Z", "Zambia", "Ethiopia")

NewMutant$population <- factor(NewMutant$population)

# removing wings from flies raised at 5% nutrition
NewMutant <- subset(NewMutant, subset = nutrition != "f5", drop = T)
NewMutant$nutrition <- droplevels(NewMutant$nutrition)
NewMutant$condition <- droplevels(NewMutant$condition)

# creating numeric and zeroed temperature variables
new_temp <- NewMutant$temperature
levels(new_temp) <- c("24", "18", "28")
temp_num <- as.numeric(levels(new_temp))[new_temp]
NewMutant <- cbind(NewMutant, temp_num)
NewMutant$temp_zeroed <- NewMutant$temp_num - 18



# Model testing the effects of temperature sex and population on frequency of wing defects
mut_freq.lmm4 <- glmer(Mutant ~ (temperature + sex + population)^3 +
                         (1 + temperature|line:population), family = binomial,
                       data = subset(NewMutant, nutrition == "f100"))
summary(mut_freq.lmm4) 

# Supplementary Data S20
Anova(mut_freq.lmm4)

# extracting marginal means
TempMut_marginal_means <- as.data.frame(Effect(c("population", "temperature" ,"sex"), mut_freq.lmm4))

# takling the average frequency of wing defects for each line for plotting
PropMutantsByLine <- aggregate(NewMutant$Mutant,
                               c(NewMutant["sex"],
                                 NewMutant["line"],
                                 NewMutant["population"],
                                 NewMutant["nutrition"],
                                 NewMutant["temperature"],
                                 NewMutant["temp_zeroed"]), 
                               mean, na.rm=TRUE )
# fixing levels
levels(PropMutantsByLine$sex) <- c("F", "M")
levels(PropMutantsByLine$line) <- c("e112", "e117", "e119", "e122", "e131", "e136", "e15",  "e16",  "e19",  "e39", "e54",  "e59",  "e65",  "e73",  "e98",  "z124", "z159", "z160", "z186", "z216","z251", "z254", "z311", "z360", "z366", "z367", "z455", "z461")

# takling the average frequency of wing defects for each population for plotting
PropMutantsByPopulation <- aggregate(NewMutant$Mutant,
                                     c(NewMutant["sex"],
                                       NewMutant["population"],
                                       NewMutant["nutrition"],
                                       NewMutant["temperature"],
                                       NewMutant["temp_zeroed"]), 
                                     mean, na.rm=TRUE )


# Supplementary Figure S9A
ag <- ggplot(data = TempMut_marginal_means, 
             aes(x = temperature, y = fit, 
                 colour = population, 
                 shape = sex)) +
  geom_point(size = 14, position = position_dodge(0.6), show.legend = FALSE) +
  geom_errorbar(aes(ymin=lower, ymax=upper), 
              lwd = 1.5, width = 0.6, position = position_dodge(0.6), show.legend = FALSE) +
  geom_jitter(data = PropMutantsByLine, aes(x = temperature, y = x), 
              size = 6, show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), alpha = 0.6)+
  theme(legend.position = c(0.70, 0.9),
        legend.box = "horizontal",
        legend.key.width = unit(1.5,"cm"),
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  scale_x_discrete(name = expression("Temperature " ( degree*C)), 
                     labels = c(18, 24, 28)) +
  ylab("Proportion of Defects")+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name = "Sex",
                       labels = c("Female", "Male"))
ag

# For legend
ag_leg <- ggplot(data = TempMut_marginal_means, 
             aes(x = temperature, y = fit, 
                 colour = population, 
                 shape = sex)) +
  geom_point(size = 14, position = position_dodge(0.6), show.legend = TRUE) +
  geom_errorbar(aes(ymin=lower, ymax=upper), 
                lwd = 1.5, width = 0.6, position = position_dodge(0.6), show.legend = FALSE) +
  geom_jitter(data = PropMutantsByLine, aes(x = temperature, y = x), 
              size = 6, show.legend = FALSE,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6))+
  theme(legend.position = c(0.70, 0.9),
        legend.box = "horizontal",
        legend.key.width = unit(1.5,"cm"),
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  scale_x_discrete(name = expression("Temperature " ( degree*C)), 
                   labels = c(18, 24, 28)) +
  ylab("Proportion of Defects")+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name = "Sex",
                       labels = c("Female", "Male"))
ag_leg

levels(PropMutantsByLine$sex) <- c("f", "m") # changing the sex levels to lower case levels
###############################################################

################# Temperatrure Wing Size CV ###################

############ Calculating CV/LevDev by removing sex effect first #############

# creating an interaction term for line and temperature
wings_temp$line_by_temp <- as.factor(with(wings_temp, interaction(line, temperature)))

wings_temp$line_by_temp <- droplevels(wings_temp$line_by_temp)

# empty vectors for CV and Lev Dev measures for loop
cv_out <- rep(NA, nlevels(wings_temp$line_by_temp))
ls1_out <- rep(NA, nlevels(wings_temp$line_by_temp)) # median of median LS on log transformed
ls2_out <- rep(NA, nlevels(wings_temp$line_by_temp)) # mean of median LS on log transformed


# loop to calculate CV and Lev Dev for each line at each temperature
for (i in 1:nlevels(wings_temp$line_by_temp)) {
  data_loop <- wings_temp[wings_temp$line_by_temp == levels(wings_temp$line_by_temp)[i], ]
  mod_loop <- lm(CS ~ 1 + sex, data = data_loop)
  line_mean <- (coef(mod_loop)[1] + sum(coef(mod_loop)))/2
  mod_adj_size <- resid(mod_loop) + line_mean
  cv_out[i] <- sd(mod_adj_size)/line_mean
  log_mod_adj_size <- log(mod_adj_size)
  ls1_out[i] <- median(abs(log_mod_adj_size - median(log_mod_adj_size)))
  ls2_out[i] <- mean(abs(log_mod_adj_size - median(log_mod_adj_size)))
  
  rm(data_loop, mod_loop, line_mean)
}

# Creating data frame with CV and Lev Dev measures for each line at each temperature
cv_dat <- data.frame(line.temp = levels(wings_temp$line_by_temp), cv_out, ls1_out, ls2_out)

# Adding a populaton factor
cv_dat$population <- ifelse(grepl("e", cv_dat$line.temp), "Ethiopia", "Zambia")
cv_dat$population <- as.factor(cv_dat$population)

# Separating line and temperature factors
cv_dat2 <- separate(data = cv_dat, col = line.temp, into= c("line", "temperature"), sep = "\\.")
cv_dat2[,1:2] <- lapply(cv_dat2[,1:2], factor)


# Modeling Sex correcrted CV
lmm.CVTemp_sexcorrected_gamma <- glmer(cv_out ~ population*temperature + (1|line:population),
                          data = cv_dat2, family = "Gamma",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))

summary(lmm.CVTemp_sexcorrected_gamma)
# Supplementary Table S16
Anova(lmm.CVTemp_sexcorrected_gamma)

lmm.CVTemp_sexcorrected_gamma_marginal_means <- as.data.frame(Effect(c("population", "temperature"),
                                                                     lmm.CVTemp_sexcorrected_gamma))

ah_sex_corrected <- ggplot(lmm.CVTemp_sexcorrected_gamma_marginal_means, 
             aes(x = temperature, y=fit, colour = population))+
  geom_point(size = 14, show.legend = FALSE, position = position_dodge(0.8))+
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), lwd = 1.5, position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = cv_dat2, aes(x = temperature, y = cv_out, colour = population), size = 6, position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "vertical", legend.title = element_text(size = 36), legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  scale_y_continuous(name="Wing Size Variation (CV)")+
  scale_x_discrete(name = expression("Temperature " ( degree*C)), 
                   labels = c(18, 24, 28))+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))
ah_sex_corrected

#Modeling sex corrected Lev Dev
lmm.LevDevTemp_sexcorrected_gamma <- glmer(ls1_out ~ population*temperature + (1|line:population),
                                       data = cv_dat2, family = "Gamma",
                                       control=glmerControl(optimizer="bobyqa",
                                                            optCtrl=list(maxfun=2e5)))

summary(lmm.LevDevTemp_sexcorrected_gamma)
# Supplementary Table S17
Anova(lmm.LevDevTemp_sexcorrected_gamma)

lmm.LevDevTemp_sexcorrected_gamma_marginal_means <- as.data.frame(Effect(c("population", "temperature"),
                                                                         lmm.LevDevTemp_sexcorrected_gamma))


ai_sex_corrected <- ggplot(lmm.LevDevTemp_sexcorrected_gamma_marginal_means, 
                           aes(x = temperature, y=fit, colour = population))+
  geom_point(size = 14, show.legend = FALSE, position = position_dodge(0.8))+
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), lwd = 1.5, position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = cv_dat2, aes(x = temperature, y = ls1_out, colour = population), size = 6, position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.8), legend.box = "vertical", legend.title = element_text(size = 36), legend.text = element_text(size = 38), legend.background = element_rect(fill="transparent"))+
  scale_y_continuous(name="Wing Size Variation (LD)")+
  scale_x_discrete(name = expression("Temperature " ( degree*C)), 
                   labels = c(18, 24, 28))+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))
ai_sex_corrected
###############################################################

################### Interlude: For internal checks only ################
#Check correlations of CV between lines raised at different temperauters
# Calculating correlations between CV at different temperatures for each population
cv_dat2_cvonly <- cv_dat2[,-c(4:5)]

cv_dat2_cvonly_wide <- spread(cv_dat2_cvonly, temperature, cv_out)

cv_dat2_cvonly_wide_Eth <- subset(cv_dat2_cvonly_wide, population == "Ethiopia")
cv_dat2_cvonly_wide_Zam <- subset(cv_dat2_cvonly_wide, population == "Zambia")


cor_plotCV_18_24 <- ggplot(data = cv_dat2_cvonly_wide, aes(x = t18, y = t24, colour = population))+
  geom_point(size = 14) +
  geom_smooth(method='lm',formula=y~x)
cor_plotCV_18_24
cor.test(cv_dat2_cvonly_wide_Eth$t18, cv_dat2_cvonly_wide_Eth$t24)
cor.test(cv_dat2_cvonly_wide_Zam$t18, cv_dat2_cvonly_wide_Zam$t24)

cor_plotCV_18_28 <- ggplot(data = cv_dat2_cvonly_wide, aes(x = t18, y = t28, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotCV_18_28

cor.test(cv_dat2_cvonly_wide_Eth$t18, cv_dat2_cvonly_wide_Eth$t28)
cor.test(cv_dat2_cvonly_wide_Zam$t18, cv_dat2_cvonly_wide_Zam$t28)


cor_plotCV_24_28 <- ggplot(data = cv_dat2_cvonly_wide, aes(x = t24, y = t28, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotCV_24_28
cor.test(cv_dat2_cvonly_wide_Eth$t24, cv_dat2_cvonly_wide_Eth$t28)
cor.test(cv_dat2_cvonly_wide_Zam$t24, cv_dat2_cvonly_wide_Zam$t28)


# Calculating correlations between Lev Dev at different temperatures for each population

cv_dat2_lsonly <- cv_dat2[,-c(3,5)]

cv_dat2_lsonly_wide <- spread(cv_dat2_lsonly, temperature, ls1_out)

cv_dat2_lsonly_wide_Eth <- subset(cv_dat2_lsonly_wide, population == "Ethiopia")
cv_dat2_lsonly_wide_Zam <- subset(cv_dat2_lsonly_wide, population == "Zambia")


cor_plotCV_18_24 <- ggplot(data = cv_dat2_lsonly_wide, aes(x = t18, y = t24, colour = population))+
  geom_point(size = 14) +
  geom_smooth(method='lm',formula=y~x)
cor_plotCV_18_24
cor.test(cv_dat2_lsonly_wide_Eth$t18, cv_dat2_lsonly_wide_Eth$t24)
cor.test(cv_dat2_lsonly_wide_Zam$t18, cv_dat2_lsonly_wide_Zam$t24)

cor_plotCV_18_28 <- ggplot(data = cv_dat2_lsonly_wide, aes(x = t18, y = t28, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotCV_18_28

cor.test(cv_dat2_lsonly_wide_Eth$t18, cv_dat2_lsonly_wide_Eth$t28)
cor.test(cv_dat2_lsonly_wide_Zam$t18, cv_dat2_lsonly_wide_Zam$t28)


cor_plotCV_24_28 <- ggplot(data = cv_dat2_lsonly_wide, aes(x = t24, y = t28, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotCV_24_28
cor.test(cv_dat2_lsonly_wide_Eth$t24, cv_dat2_lsonly_wide_Eth$t28)
cor.test(cv_dat2_lsonly_wide_Zam$t24, cv_dat2_lsonly_wide_Zam$t28)
################################################


################## Wingshape Total variance, Eccentricity, rSDE, rSDE2 ###################
wingstempPC$nutrition <- droplevels(wingstempPC$nutrition)
wingstempPC$condition <- droplevels(wingstempPC$condition)
wingstempPC$line <- droplevels(wingstempPC$line)

#each temperature seperately:
wingstempPC_18deg <- subset(wingstempPC, temperature == "t18", drop = TRUE)
wingstempPC_18deg$temperature <- droplevels(wingstempPC_18deg$temperature)
wingstempPC_18deg$line <- droplevels(wingstempPC_18deg$line)

wingstempPC_24deg <- subset(wingstempPC, temperature == "t24", drop = TRUE)
wingstempPC_24deg$temperature <- droplevels(wingstempPC_24deg$temperature)
wingstempPC_24deg$line <- droplevels(wingstempPC_24deg$line)

wingstempPC_28deg <- subset(wingstempPC, temperature == "t28", drop = TRUE)
wingstempPC_28deg$temperature <- droplevels(wingstempPC_28deg$temperature)
wingstempPC_28deg$line <- droplevels(wingstempPC_28deg$line)

# Calculating TV, Ecc, rSDE and rSDE for flies raised at 18 deg:

# From Annat's code (Sept 2016)
Xall.res.temp.18deg <- matrix(NA, nrow = nrow(wingstempPC_18deg), ncol = 58) # ncol is the number of PCs we want to use
line_values.temp.18deg <- rep(NA, length = nrow(wingstempPC_18deg)) #  A vector with line names.

# loging centroid size, removing variation due to sex and log centroid size, extracring residuals for within-line measures of varation calculations below. 
for (ln in levels(wingstempPC_18deg$line)) {
  ind <- which(wingstempPC_18deg$line==ln)
  Xln <- as.matrix(wingstempPC_18deg[ind,115:172]) 
  lcs <- log(wingstempPC_18deg$CS[ind]) # log centroid size
  sex <- wingstempPC_18deg$sex[ind]
  reg <- lm(Xln ~ sex + lcs) # regression model
  Xm <- colMeans(predict(reg)) #CHECK IF ITS PREDICTING THE RIGHT WAY # predicted mean configuration at mean centroid size to add back to the residuals. Predict is fitting it by default at the mean for lcs.
  Xall.res.temp.18deg[ind,] <- reg$residuals + matrix(Xm, nr=length(ind), nc=length(Xm), byrow=TRUE) # recentered on the line mean
  line_values.temp.18deg[ind] = ln
}

# centering shape
centered.shape.temp_18deg <- data.frame(wingstempPC_18deg$line, Xall.res.temp.18deg)

# Calculating TV, Ecc, rSDE and rSDE2
tot_var_temp_18deg <- rep(NA, nlevels(centered.shape.temp_18deg$wingstempPC_18deg.line))
eccentricity_temp_18deg <- rep(NA, nlevels(centered.shape.temp_18deg$wingstempPC_18deg.line))
rSDE_temp_18deg <- rep(NA, nlevels(centered.shape.temp_18deg$wingstempPC_18deg.line))
rSDE2_temp_18deg <- rep(NA, nlevels(centered.shape.temp_18deg$wingstempPC_18deg.line))

for (ln in 1:nlevels(centered.shape.temp_18deg$wingstempPC_18deg.line)) {
  lev <- levels(centered.shape.temp_18deg$wingstempPC_18deg.line)[ln]
  mat = centered.shape.temp_18deg[centered.shape.temp_18deg$wingstempPC_18deg.line == lev, 2:59]
  p <- ncol(mat) # number of variables
  cov_mat <- cov(mat)
  eig_mat <- eigen(cov_mat)$values
  tot_var_temp_18deg[ln] <- sum(eig_mat) 
  
  eccentricity_temp_18deg[ln] <- eig_mat[1]/(sum(eig_mat))
  
  rSDE_temp_18deg[ln] <- sqrt(var(eig_mat)*((p-1)/p)) # rescaled to p instead of p - 1 
  
  
  rSDE2_temp_18deg[ln] <- sqrt(var(eig_mat/sum(eig_mat))*((p-1)/p) )# rescaling this by the total variance. I think since we have already applied some corrections for size and sex this may result in over-correction.
}

# putting all the measures into a data frame
integration_measures_temp_18deg <- data.frame(tot_var_temp_18deg, 
                                              eccentricity_temp_18deg, 
                                              rSDE_temp_18deg, 
                                              rSDE2_temp_18deg, 
                                              line = levels(centered.shape.temp_18deg$wingstempPC_18deg.line))

# creating a population factor for data frame
integration_measures_temp_18deg$pop <- as.factor(c(rep("Ethiopia", 10), rep("Zambia", 9)))



# Calculating TV, Ecc, rSDE and rSDE for flies raised at 24 deg:
# From Annat's code (Sept 2016)
Xall.res.temp.24deg <- matrix(NA, nrow = nrow(wingstempPC_24deg), ncol = 58) # ncol is the number of PCs we want to use
line_values.temp.24deg <- rep(NA, length = nrow(wingstempPC_24deg)) #  A vector with line names.

# loging centroid size, removing variation due to sex and log centroid size, extracring residuals for within-line measures of varation calculations below.
for (ln in levels(wingstempPC_24deg$line)) {
  ind <- which(wingstempPC_24deg$line==ln)
  Xln <- as.matrix(wingstempPC_24deg[ind,115:172]) 
  lcs <- log(wingstempPC_24deg$CS[ind]) # log centroid size
  sex <- wingstempPC_24deg$sex[ind]
  reg <- lm(Xln ~ sex + lcs) # regression model
  Xm <- colMeans(predict(reg)) #CHECK IF ITS PREDICTING THE RIGHT WAY # predicted mean configuration at mean centroid size to add back to the residuals. Predict is fitting it by default at the mean for lcs.
  Xall.res.temp.24deg[ind,] <- reg$residuals + matrix(Xm, nr=length(ind), nc=length(Xm), byrow=TRUE) # recentered on the line mean
  line_values.temp.24deg[ind] = ln
}

# centering shape data
centered.shape.temp_24deg <- data.frame(wingstempPC_24deg$line, Xall.res.temp.24deg)

# Calculating TV, Ecc, rSDE and rSDE2
tot_var_temp_24deg <- rep(NA, nlevels(centered.shape.temp_24deg$wingstempPC_24deg.line))
eccentricity_temp_24deg <- rep(NA, nlevels(centered.shape.temp_24deg$wingstempPC_24deg.line))
rSDE_temp_24deg <- rep(NA, nlevels(centered.shape.temp_24deg$wingstempPC_24deg.line))
rSDE2_temp_24deg <- rep(NA, nlevels(centered.shape.temp_24deg$wingstempPC_24deg.line))


for (ln in 1:nlevels(centered.shape.temp_24deg$wingstempPC_24deg.line)) {
  lev <- levels(centered.shape.temp_24deg$wingstempPC_24deg.line)[ln]
  mat = centered.shape.temp_24deg[centered.shape.temp_24deg$wingstempPC_24deg.line == lev, 2:59]
  p <- ncol(mat) # number of variables
  cov_mat <- cov(mat)
  eig_mat <- eigen(cov_mat)$values
  tot_var_temp_24deg[ln] <- sum(eig_mat) 
  
  eccentricity_temp_24deg[ln] <- eig_mat[1]/(sum(eig_mat))
  
  rSDE_temp_24deg[ln] <- sqrt(var(eig_mat)*((p-1)/p)) # rescaled to p instead of p - 1 
  
  
  rSDE2_temp_24deg[ln] <- sqrt(var(eig_mat/sum(eig_mat))*((p-1)/p) )# rescaling this by the total variance. I think since we have already applied some corrections for size and sex this may result in over-correction.
}

# putting all the measures into a data frame
integration_measures_temp_24deg <- data.frame(tot_var_temp_24deg, 
                                              eccentricity_temp_24deg, 
                                              rSDE_temp_24deg, 
                                              rSDE2_temp_24deg, 
                                              line = levels(centered.shape.temp_24deg$wingstempPC_24deg.line))

# creating a population factor for data frame
integration_measures_temp_24deg$pop <- as.factor(c(rep("Ethiopia", 13), rep("Zambia", 9)))


# Calculating TV, Ecc, rSDE and rSDE for flies raised at 28 deg:

# From Annat's code (Sept 2016)
Xall.res.temp.28deg <- matrix(NA, nrow = nrow(wingstempPC_28deg), ncol = 58) # ncol is the number of PCs we want to use
line_values.temp.28deg <- rep(NA, length = nrow(wingstempPC_28deg)) #  A vector with line names.

# loging centroid size, removing variation due to sex and log centroid size, extracring residuals for within-line measures of varation calculations below.
for (ln in levels(wingstempPC_28deg$line)) {
  ind <- which(wingstempPC_28deg$line==ln)
  Xln <- as.matrix(wingstempPC_28deg[ind,115:172]) 
  lcs <- log(wingstempPC_28deg$CS[ind]) # log centroid size
  sex <- wingstempPC_28deg$sex[ind]
  reg <- lm(Xln ~ sex + lcs) # regression model
  Xm <- colMeans(predict(reg)) #CHECK IF ITS PREDICTING THE RIGHT WAY # predicted mean configuration at mean centroid size to add back to the residuals. Predict is fitting it by default at the mean for lcs.
  Xall.res.temp.28deg[ind,] <- reg$residuals + matrix(Xm, nr=length(ind), nc=length(Xm), byrow=TRUE) # recentered on the line mean
  line_values.temp.28deg[ind] = ln
}

#centering shape data
centered.shape.temp_28deg <- data.frame(wingstempPC_28deg$line, Xall.res.temp.28deg)

# Calculating TV, Ecc, rSDE and rSDE2
tot_var_temp_28deg <- rep(NA, nlevels(centered.shape.temp_28deg$wingstempPC_28deg.line))
eccentricity_temp_28deg <- rep(NA, nlevels(centered.shape.temp_28deg$wingstempPC_28deg.line))
rSDE_temp_28deg <- rep(NA, nlevels(centered.shape.temp_28deg$wingstempPC_28deg.line))
rSDE2_temp_28deg <- rep(NA, nlevels(centered.shape.temp_28deg$wingstempPC_28deg.line))


for (ln in 1:nlevels(centered.shape.temp_28deg$wingstempPC_28deg.line)) {
  lev <- levels(centered.shape.temp_28deg$wingstempPC_28deg.line)[ln]
  mat = centered.shape.temp_28deg[centered.shape.temp_28deg$wingstempPC_28deg.line == lev, 2:59]
  p <- ncol(mat) # number of variables
  cov_mat <- cov(mat)
  eig_mat <- eigen(cov_mat)$values
  tot_var_temp_28deg[ln] <- sum(eig_mat) 
  
  eccentricity_temp_28deg[ln] <- eig_mat[1]/(sum(eig_mat))
  
  rSDE_temp_28deg[ln] <- sqrt(var(eig_mat)*((p-1)/p)) # rescaled to p instead of p - 1 
  
  
  rSDE2_temp_28deg[ln] <- sqrt(var(eig_mat/sum(eig_mat))*((p-1)/p) )# rescaling this by the total variance. I think since we have already applied some corrections for size and sex this may result in over-correction.
}

# putting all the measures into a data frame
integration_measures_temp_28deg <- data.frame(tot_var_temp_28deg, 
                                              eccentricity_temp_28deg, 
                                              rSDE_temp_28deg, 
                                              rSDE2_temp_28deg, 
                                              line = levels(centered.shape.temp_28deg$wingstempPC_28deg.line))

# creating a population factor for data frame
integration_measures_temp_28deg$pop <- as.factor(c(rep("Ethiopia", 10), rep("Zambia", 8)))



# Fixing the data frames for the within-line measures of variation for each temperature so that we can consolidate them into a single data frame.
names(integration_measures_temp_18deg) <- c("tot_var_temp", 
                                            "eccentricity_temp", 
                                            "rSDE_temp",
                                            "rSDE2_temp",
                                            "line",
                                            "pop")
integration_measures_temp_18deg$temperature <- as.factor(c(rep("t18",19)))

names(integration_measures_temp_24deg) <- c("tot_var_temp", 
                                            "eccentricity_temp", 
                                            "rSDE_temp",
                                            "rSDE2_temp",
                                            "line",
                                            "pop")
integration_measures_temp_24deg$temperature <- as.factor(c(rep("t24",22)))

names(integration_measures_temp_28deg) <- c("tot_var_temp", 
                                            "eccentricity_temp", 
                                            "rSDE_temp",
                                            "rSDE2_temp",
                                            "line",
                                            "pop")
integration_measures_temp_28deg$temperature <- as.factor(c(rep("t28",18)))

# consolidating the data frames of within-line measures of variation for shape into a single data frame
integration_measures_temp <- rbind(integration_measures_temp_18deg,
                                   integration_measures_temp_24deg,
                                   integration_measures_temp_28deg)

# removing frequency of mutational defects data from flies raised at 15% nutrition
PropMutantsbyLineTemp <- subset(PropMutantsByLine, nutrition == "f100")

# fixing column names of frequency of mutational defects data frame for comparisons later
names(PropMutantsbyLineTemp) <- c("sex","line", "population", "nutrition", "temperature", "temp_zeroed", "prop_mut")

# fixing column names of within-line variation for wing shape data frame
names(integration_measures_temp) <- c("tot_var_temp", "eccentricity_temp", "rSDE_temp","rSDE2_temp", "line", "population", "temperature")


#glmer model with gamma distribution

#total variance multiplied by a factor of 1000

integration_measures_temp$tot_var_temp_alt <- integration_measures_temp$tot_var_temp*1000

# Generalized linear mixed model for Total variance
lmm.TVTemp_gamma <- glmer(tot_var_temp_alt ~ population*temperature + (1|line:population),
                          data = integration_measures_temp, family = "Gamma",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))

summary(lmm.TVTemp_gamma)


#Supplementary Table S18
Anova(lmm.TVTemp_gamma)

lmm.TVTemp_gamma_marginal_means <- as.data.frame(Effect(c("population", "temperature"),
                                                        lmm.TVTemp_gamma))

# Figure 5C
aj <- ggplot(lmm.TVTemp_gamma_marginal_means, 
              aes(x = temperature, y=fit, colour = population))+
  geom_point(size = 14, show.legend = FALSE, position = position_dodge(0.8))+
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), lwd = 1.5, position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = integration_measures_temp, aes(x = temperature, y = tot_var_temp_alt), size = 6, position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.9), legend.box = "horizontal")+
  scale_y_continuous(name="Wing Shape Variation (TV)")+
  scale_x_discrete(name = expression("Temperature " ( degree*C)), 
                   labels = c(18, 24, 28))+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))
aj

# for extractring legend
aj_leg<- ggplot(lmm.TVTemp_gamma_marginal_means, 
              aes(x = temperature, y=fit, colour = population))+
  geom_point(size = 14, show.legend = TRUE, position = position_dodge(0.8))+
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), lwd = 1.5, position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = integration_measures_temp, aes(x = temperature, y = tot_var_temp_alt), size = 6, position = position_jitterdodge(jitter.width = 0.2), show.legend = TRUE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.9), legend.box = "horizontal")+
  scale_y_continuous(name="Wing Shape Variation (TV)")+
  scale_x_discrete(name = expression("Temperature " ( degree*C)), 
                   labels = c(18, 24, 28))+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))+
  scale_shape_discrete(name = "Sex",
                      labels = c("Female", "Male"))
aj_leg

# Generalized linear mixed model for eccentricity
lmm.EccTemp_gamma <- glmer(eccentricity_temp ~ population*temperature + (1|line:population),
                          data = integration_measures_temp, family = "Gamma",
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))

summary(lmm.EccTemp_gamma)

#Supplementary Table S19
Anova(lmm.EccTemp_gamma)

lmm.EccTemp_gamma_marginal_means <- as.data.frame(Effect(c("population", "temperature"),
                                                         lmm.EccTemp_gamma))

# Figure 5D
ak <- ggplot(lmm.EccTemp_gamma_marginal_means, 
              aes(x = temperature, y=fit, colour = population))+
  geom_point(size = 14, show.legend = FALSE, position = position_dodge(0.8))+
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), lwd = 1.5, position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = integration_measures_temp, aes(x = temperature, y = eccentricity_temp), size = 6, position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.9), legend.box = "horizontal")+
  scale_y_continuous(name="Wing Shape Integration (Ecc)")+
  scale_x_discrete(name = expression("Temperature " ( degree*C)), 
                   labels = c(18, 24, 28))+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))
ak


#rSDE multiplied by a factor of 10000
integration_measures_temp$rSDE_temp_alt <- integration_measures_temp$rSDE_temp * 10000

# Generalized linear mixed model for rSDE
lmm.rSDETemp_gamma <- glmer(rSDE_temp_alt ~ population*temperature + (1|line:population),
                           data = integration_measures_temp, family = "Gamma",
                           control=glmerControl(optimizer="bobyqa",
                                                optCtrl=list(maxfun=2e5)))

summary(lmm.rSDETemp_gamma)
Anova(lmm.rSDETemp_gamma)

lmm.rSDETemp_gamma_marginal_means <- as.data.frame(Effect(c("population", "temperature"),
                                                          lmm.rSDETemp_gamma))
# Supplementary Figure S9C
al<- ggplot(lmm.rSDETemp_gamma_marginal_means, 
              aes(x = temperature, y=fit, colour = population))+
  geom_point(size = 14, show.legend = FALSE, position = position_dodge(0.8))+
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), lwd = 1.5, position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = integration_measures_temp, aes(x = temperature, y = rSDE_temp_alt), size = 6, position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.9), legend.box = "horizontal")+
  scale_y_continuous(name="Wing Shape Integration (rSDE)")+
  scale_x_discrete(name = expression("Temperature " ( degree*C)), 
                   labels = c(18, 24, 28))+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))
al



# rSDE2 multipled by 10
integration_measures_temp$rSDE2_temp_alt <- integration_measures_temp$rSDE2_temp * 10

# Generalized linear mixed model for rSDE2
lmm.rSDE2_Temp_gamma <- glmer(rSDE2_temp_alt ~ population*temperature + (1|line:population),
                            data = integration_measures_temp, family = "Gamma",
                            control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))

summary(lmm.rSDE2_Temp_gamma)
Anova(lmm.rSDE2_Temp_gamma)

lmm.rSDE2_Temp_gamma_marginal_means <- as.data.frame(Effect(c("population", "temperature"),
                                                          lmm.rSDE2_Temp_gamma))

# Supplementary Figure S9D
am <- ggplot(lmm.rSDE2_Temp_gamma_marginal_means, 
              aes(x = temperature, y=fit, colour = population))+
  geom_point(size = 14, show.legend = FALSE, position = position_dodge(0.8))+
  geom_errorbar(width=0.3, aes(ymin=lower, ymax=upper), lwd = 1.5, position = position_dodge(0.8), show.legend = FALSE)+
  geom_jitter(data = integration_measures_temp, aes(x = temperature, y = rSDE2_temp_alt), size = 6, position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.9), legend.box = "horizontal")+
  scale_y_continuous(name="Wing Shape Integration (rSDE2)")+
  scale_x_discrete(name = expression("Temperature " ( degree*C)), 
                   labels = c(18, 24, 28))+
  scale_colour_manual(name = "Population",
                      values = c("black","gray60"),
                      labels = c("High-Altitude", "Low-Altitude"))
am

################################################################################

################ Interlude: For checking only ####################
# Looking at the correlations between within line measures of variation for wing shape for lines raised at different temperatures

#total variance
integration_measures_temp_totvar <- integration_measures_temp[,c(1,5:7)]

integration_measures_temp_totvar_wide <- spread(integration_measures_temp_totvar, temperature,  tot_var_temp)

integration_measures_temp_totvar_wide_Eth <- subset(integration_measures_temp_totvar_wide, 
                                                    population == "Ethiopia")
integration_measures_temp_totvar_wide_Zam <- subset(integration_measures_temp_totvar_wide, 
                                                    population == "Zambia")

cor_plotTV_18_24 <- ggplot(data = integration_measures_temp_totvar_wide, 
                           aes(x = t18, y = t24, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotTV_18_24
cor.test(integration_measures_temp_totvar_wide_Eth$t18, integration_measures_temp_totvar_wide_Eth$t24)
cor.test(integration_measures_temp_totvar_wide_Zam$t18, integration_measures_temp_totvar_wide_Zam$t24)


cor_plotTV_18_28 <- ggplot(data = integration_measures_temp_totvar_wide, 
                           aes(x = t18, y = t28, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotTV_18_28

cor.test(integration_measures_temp_totvar_wide_Eth$t18, integration_measures_temp_totvar_wide_Eth$t28)
cor.test(integration_measures_temp_totvar_wide_Zam$t18, integration_measures_temp_totvar_wide_Zam$t28)



cor_plotTV_18_28 <- ggplot(data = integration_measures_temp_totvar_wide, 
                           aes(x = t24, y = t28, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotTV_18_28

cor.test(integration_measures_temp_totvar_wide_Eth$t24, integration_measures_temp_totvar_wide_Eth$t28)
cor.test(integration_measures_temp_totvar_wide_Zam$t24, integration_measures_temp_totvar_wide_Zam$t28)

#eccentricity
integration_measures_temp_eccentricity <- integration_measures_temp[,c(2, 5,6,7)]

integration_measures_temp_eccentricity_wide <- spread(integration_measures_temp_eccentricity, temperature,  eccentricity_temp)

integration_measures_temp_eccentricity_wide_Eth <- subset(integration_measures_temp_eccentricity_wide,
                                                          population == "Ethiopia")

integration_measures_temp_eccentricity_wide_Zam <- subset(integration_measures_temp_eccentricity_wide,
                                                          population == "Zambia")

cor_plotEcc_18_24 <- ggplot(data = integration_measures_temp_eccentricity_wide, 
                            aes(x = t18, y = t24, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotEcc_18_24

cor.test(integration_measures_temp_eccentricity_wide_Eth$t18, integration_measures_temp_eccentricity_wide_Eth$t24)
cor.test(integration_measures_temp_eccentricity_wide_Zam$t18, integration_measures_temp_eccentricity_wide_Zam$t24)


cor_plotEcc_18_28 <- ggplot(data = integration_measures_temp_eccentricity_wide, 
                            aes(x = t18, y = t28, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotEcc_18_28

cor.test(integration_measures_temp_eccentricity_wide_Eth$t18, integration_measures_temp_eccentricity_wide_Eth$t28)
cor.test(integration_measures_temp_eccentricity_wide_Zam$t18, integration_measures_temp_eccentricity_wide_Zam$t28)


cor_plotEcc_24_28 <- ggplot(data = integration_measures_temp_eccentricity_wide, 
                            aes(x = t24, y = t28, colour = population))+
  geom_point(size = 14)+
  geom_smooth(method='lm',formula=y~x)
cor_plotEcc_24_28

cor.test(integration_measures_temp_eccentricity_wide_Eth$t24, integration_measures_temp_eccentricity_wide_Eth$t28)
cor.test(integration_measures_temp_eccentricity_wide_Zam$t24, integration_measures_temp_eccentricity_wide_Zam$t28)

##################################################################

######## Measures of Variablity temp data vs mutations ########

# Wing size CV vs Defects
Mutants_CV_lines <- merge(PropMutantsbyLineTemp, cv_dat2, by = c("population", "line", "temperature"))


# Supplementary Figure S9A
an <- ggplot(Mutants_CV_lines, aes( y = prop_mut, x = cv_out, 
                                 color = population, shape = temperature)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.8, 0.7), 
        legend.box = "vertical",
        legend.background = element_rect(fill = alpha("white", 0), size=0.5),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))+
  labs(x= "Wing Size Variation (CV)", y="Proportion of Defects")+
  scale_color_manual(name= "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  =expression("Temperature " ( degree*C)), 
                       labels = c("18 deg", "24 deg", "28 deg"))+
  scale_x_continuous(limits = c(0, 0.105), breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.1))
an


an_leg <- ggplot(Mutants_CV_lines, aes( y = prop_mut, x = cv_out, 
                                    color = population, shape = temperature)) +
  geom_point(size = 10, show.legend = TRUE)+
  theme(legend.position = c(0.8, 0.7), 
        legend.box = "vertical",
        legend.background = element_rect(fill = alpha("white", 0), size=0.5),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))+
  labs(x= "Wing Size Variation (CV)", y="Proportion of Defects")+
  scale_color_manual(name= "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  = expression("Temperature " ( degree*C)), 
                       labels = c("18 deg", "24 deg", "28 deg"))+
  scale_x_continuous(limits = c(0, 0.105), breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.1))
an_leg

Mutants_CV_lines_Eth_18 <- subset(Mutants_CV_lines, population == "Ethiopia" & temperature == "t18")
Mutants_CV_lines_Zam_18 <- subset(Mutants_CV_lines, population == "Zambia" & temperature == "t18")

Mutants_CV_lines_Eth_24 <- subset(Mutants_CV_lines, population == "Ethiopia" & temperature == "t24")
Mutants_CV_lines_Zam_24 <- subset(Mutants_CV_lines, population == "Zambia" & temperature == "t24")

Mutants_CV_lines_Eth_28 <- subset(Mutants_CV_lines, population == "Ethiopia" & temperature == "t28")
Mutants_CV_lines_Zam_28 <- subset(Mutants_CV_lines, population == "Zambia" & temperature == "t28")

cor.test(Mutants_CV_lines_Eth_18$prop_mut, Mutants_CV_lines_Eth_18$cv_out)

cor.test(Mutants_CV_lines_Zam_18$prop_mut, Mutants_CV_lines_Zam_18$cv_out)


cor.test(Mutants_CV_lines_Eth_24$prop_mut, Mutants_CV_lines_Eth_24$cv_out)

cor.test(Mutants_CV_lines_Zam_24$prop_mut, Mutants_CV_lines_Zam_24$cv_out)


cor.test(Mutants_CV_lines_Eth_28$prop_mut, Mutants_CV_lines_Eth_28$cv_out)

cor.test(Mutants_CV_lines_Zam_28$prop_mut, Mutants_CV_lines_Zam_28$cv_out)


# Wing size Levene Deviates vs proportion of defects

ao <- ggplot(Mutants_CV_lines, aes( y = prop_mut, x = ls1_out, 
                                    color = population, shape = temperature)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.7, 0.8), 
        legend.box = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0), size=0.5),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))+
  labs(x= "Wing Size Variation (LD)", y="Proportion of Defects")+
  scale_color_manual(name= "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  ="Temperature", 
                       labels = c("18 deg", "24 deg", "28 deg"))
ao


cor.test(Mutants_CV_lines$prop_mut, Mutants_CV_lines$ls1_out)

cor.test(Mutants_CV_lines$prop_mut, Mutants_CV_lines$ls1_out)



cor.test(Mutants_CV_lines_Eth_18$prop_mut, Mutants_CV_lines_Eth_18$ls1_out)

cor.test(Mutants_CV_lines_Zam_18$prop_mut, Mutants_CV_lines_Zam_18$ls1_out)


cor.test(Mutants_CV_lines_Eth_24$prop_mut, Mutants_CV_lines_Eth_24$ls1_out)

cor.test(Mutants_CV_lines_Zam_24$prop_mut, Mutants_CV_lines_Zam_24$ls1_out)


cor.test(Mutants_CV_lines_Eth_28$prop_mut, Mutants_CV_lines_Eth_28$ls1_out)

cor.test(Mutants_CV_lines_Zam_28$prop_mut, Mutants_CV_lines_Zam_28$ls1_out)





PropMutantsbyLineTemp_averagedMF <- aggregate(PropMutantsbyLineTemp["prop_mut"],
                                              by = c(PropMutantsbyLineTemp["line"], 
                                                     PropMutantsbyLineTemp["population"],
                                                     PropMutantsbyLineTemp["temperature"]),
                                              mean, na.rm = TRUE)


# Wing shape Total variance and Eccentricity vs proportion of defects
Mutants_Integ_lines <- merge(PropMutantsbyLineTemp_averagedMF, integration_measures_temp, by = c("population", "line", "temperature"))

# Supplementary Figure S10B
ap <- ggplot(Mutants_Integ_lines, aes( y = prop_mut, x = tot_var_temp_alt, 
                                        color = population, shape = temperature)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.7, 0.8), 
        legend.box = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0), size=0.5),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))+
  labs(x= "Wing Shape Variation (TV)", y="Proportion of Defects")+
  scale_color_manual(name= "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  ="Temperature", 
                        labels = c("18 deg", "24 deg", "28 deg"))+
  scale_x_continuous(limits = c(0.08, 0.31), breaks = c( 0.10, 0.15, 0.20, 0.25, 0.30))
ap


Mutants_Integ_lines_Eth_18 <- subset(Mutants_Integ_lines, population == "Ethiopia" & temperature == "t18")
Mutants_Integ_lines_Zam_18 <- subset(Mutants_Integ_lines, population == "Zambia" & temperature == "t18")

Mutants_Integ_lines_Eth_24 <- subset(Mutants_Integ_lines, population == "Ethiopia" & temperature == "t24")
Mutants_Integ_lines_Zam_24 <- subset(Mutants_Integ_lines, population == "Zambia" & temperature == "t24")

Mutants_Integ_lines_Eth_28 <- subset(Mutants_Integ_lines, population == "Ethiopia" & temperature == "t28")
Mutants_Integ_lines_Zam_28 <- subset(Mutants_Integ_lines, population == "Zambia" & temperature == "t28")

cor.test(Mutants_Integ_lines_Eth_18$prop_mut, Mutants_Integ_lines_Eth_18$tot_var_temp_alt)
cor.test(Mutants_Integ_lines_Zam_18$prop_mut, Mutants_Integ_lines_Zam_18$tot_var_temp_alt)

cor.test(Mutants_Integ_lines_Eth_24$prop_mut, Mutants_Integ_lines_Eth_24$tot_var_temp_alt)
cor.test(Mutants_Integ_lines_Zam_24$prop_mut, Mutants_Integ_lines_Zam_24$tot_var_temp_alt)

cor.test(Mutants_Integ_lines_Eth_28$prop_mut, Mutants_Integ_lines_Eth_28$tot_var_temp_alt)
cor.test(Mutants_Integ_lines_Zam_28$prop_mut, Mutants_Integ_lines_Zam_28$tot_var_temp_alt)


#Supplementary Figure S10C
aq <- ggplot(Mutants_Integ_lines, aes( y = prop_mut, x = eccentricity_temp, 
                                       color = population, shape = temperature)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.7, 0.8), 
        legend.box = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0), size=0.5),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))+
  labs(x= "Wing Shape Integration (Ecc)", y="Proportion of Defects")+
  scale_color_manual(name= "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  ="Temperature", 
                       labels = c("18 deg", "24 deg", "28 deg"))
aq


cor.test(Mutants_Integ_lines_Eth_18$prop_mut, Mutants_Integ_lines_Eth_18$eccentricity_temp)
cor.test(Mutants_Integ_lines_Zam_18$prop_mut, Mutants_Integ_lines_Zam_18$eccentricity_temp)

cor.test(Mutants_Integ_lines_Eth_24$prop_mut, Mutants_Integ_lines_Eth_24$eccentricity_temp)
cor.test(Mutants_Integ_lines_Zam_24$prop_mut, Mutants_Integ_lines_Zam_24$eccentricity_temp)

cor.test(Mutants_Integ_lines_Eth_28$prop_mut, Mutants_Integ_lines_Eth_28$eccentricity_temp)
cor.test(Mutants_Integ_lines_Zam_28$prop_mut, Mutants_Integ_lines_Zam_28$eccentricity_temp)


#Supplementary Figure S10D
ar <- ggplot(Mutants_Integ_lines, aes( y = prop_mut, x = rSDE_temp_alt, 
                                       color = population, shape = temperature)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.7, 0.8), 
        legend.box = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0), size=0.5),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))+
  labs(x= "Wing Shape Integration (rSDE)", y="Proportion of Defects")+
  scale_color_manual(name= "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  ="Temperature", 
                       labels = c("18 deg", "24 deg", "28 deg"))
ar


cor.test(Mutants_Integ_lines_Eth_18$prop_mut, Mutants_Integ_lines_Eth_18$rSDE_temp_alt)
cor.test(Mutants_Integ_lines_Zam_18$prop_mut, Mutants_Integ_lines_Zam_18$rSDE_temp_alt)

cor.test(Mutants_Integ_lines_Eth_24$prop_mut, Mutants_Integ_lines_Eth_24$rSDE_temp_alt)
cor.test(Mutants_Integ_lines_Zam_24$prop_mut, Mutants_Integ_lines_Zam_24$rSDE_temp_alt)

cor.test(Mutants_Integ_lines_Eth_28$prop_mut, Mutants_Integ_lines_Eth_28$rSDE_temp_alt)
cor.test(Mutants_Integ_lines_Zam_28$prop_mut, Mutants_Integ_lines_Zam_28$rSDE_temp_alt)


#rSDE2 vs defects
ar_2 <- ggplot(Mutants_Integ_lines, aes( y = prop_mut, x = rSDE2_temp, 
                                       color = population, shape = temperature)) +
  geom_point(size = 10, show.legend = FALSE)+
  theme(legend.position = c(0.7, 0.8), 
        legend.box = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0), size=0.5),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))+
  labs(x= "Wing Shape Integration (rSDE2)", y="Proportion of Defects")+
  scale_color_manual(name= "Population", 
                     labels = c("High-Altitude","Low-Altitude"), 
                     values=c("black","gray60"))+
  scale_shape_discrete(name  ="Temperature", 
                       labels = c("18 deg", "24 deg", "28 deg"))
ar_2


cor.test(Mutants_Integ_lines_Eth_18$prop_mut, Mutants_Integ_lines_Eth_18$rSDE2_temp)
cor.test(Mutants_Integ_lines_Zam_18$prop_mut, Mutants_Integ_lines_Zam_18$rSDE2_temp)

cor.test(Mutants_Integ_lines_Eth_24$prop_mut, Mutants_Integ_lines_Eth_24$rSDE2_temp)
cor.test(Mutants_Integ_lines_Zam_24$prop_mut, Mutants_Integ_lines_Zam_24$rSDE2_temp)

cor.test(Mutants_Integ_lines_Eth_28$prop_mut, Mutants_Integ_lines_Eth_28$rSDE2_temp)
cor.test(Mutants_Integ_lines_Zam_28$prop_mut, Mutants_Integ_lines_Zam_28$rSDE2_temp)

###############################################################


###################### Wing Shape Plots ########################

# producing mean shapes for each pop and temp
temp_wings_shape_means<- aggregate( wings_temp[,7:103], 
                                    by=list( sex=wings_temp$sex,
                                             population=wings_temp$population,
                                             temperature = wings_temp$temperature),
                                    FUN=mean)



procoords <- temp_wings_shape_means[,4:99] * matrix( rep( rep( c(-1, 1), 48), 12), 
                                        ncol= 96, byrow = T)



norm_vec <- function(x) sqrt(sum(x^2))


# Ethiopia vs Zambia at 18 deg Males
Eth_Zam_18_Males_dif <-  procoords[4, ] -  procoords[2, ]
Eth_Zam_18_Males_mean <- colMeans(procoords[c(4,2),])
EthZam_18_Males_PD <- norm_vec(Eth_Zam_18_Males_dif)

pdf(file="../outputs/WingEffectPlots/wingeffect_Eth_Zam_t18_males.pdf")
wingEffect_EthZam_t18_Males <- WingEffect( Eth_Zam_18_Males_mean, Eth_Zam_18_Males_dif , 
                                     Eth_Zam_18_Males_dif,
                                     wingcol=c("black","black", "gray60" ),
                                     scale.factor = 2,
                                     scale.display = FALSE,
                                     wingframe = FALSE,
                                     winglwd=c(3, 3, 3))
dev.off()


Eth_Zam_18_Males_dif

# Ethiopia vs Zambia at 24 deg Males
Eth_Zam_24_Males_dif <-  procoords[8, ] -  procoords[6, ]
Eth_Zam_24_Males_mean <- colMeans(procoords[c(8,6),])
EthZam_24_Males_PD <- norm_vec(Eth_Zam_24_Males_dif)

pdf(file="../outputs/WingEffectPlots/wingeffect_Eth_Zam_t24_males.pdf")
wingEffect_EthZam_t24_Males <- WingEffect( Eth_Zam_24_Males_mean, Eth_Zam_24_Males_dif , 
                                           Eth_Zam_24_Males_dif,
                                           wingcol=c("black","black", "gray60" ),
                                           scale.factor = 2,
                                           scale.display = FALSE,
                                           wingframe = FALSE,
                                           winglwd=c(3, 3, 3))
dev.off()

# Ethiopia vs Zambia at 28 deg Males
Eth_Zam_28_Males_dif <-  procoords[12, ] -  procoords[10, ]
Eth_Zam_28_Males_mean <- colMeans(procoords[c(12,10),])
EthZam_28_Males_PD <- norm_vec(Eth_Zam_28_Males_dif)

pdf(file="../outputs/WingEffectPlots/wingeffect_Eth_Zam_t28_males.pdf")
wingEffect_EthZam_t24_Males <- WingEffect( Eth_Zam_28_Males_mean, Eth_Zam_28_Males_dif , 
                                           Eth_Zam_28_Males_dif,
                                           wingcol=c("black","black", "gray60" ),
                                           scale.factor = 2,
                                           scale.display = FALSE,
                                           wingframe = FALSE,
                                           winglwd=c(3, 3, 3))
dev.off()




# Ethiopia vs Zambia at 18 deg Females
Eth_Zam_18_Females_dif <-  procoords[3, ] -  procoords[1, ]
Eth_Zam_18_Females_mean <- colMeans(procoords[c(3,1),])
EthZam_18_Females_PD <- norm_vec(Eth_Zam_18_Females_dif)

pdf(file="../outputs/WingEffectPlots/wingeffect_Eth_Zam_t18_Females.pdf")
wingEffect_EthZam_t18_Females <- WingEffect( Eth_Zam_18_Females_mean, Eth_Zam_18_Females_dif , 
                                           Eth_Zam_18_Females_dif,
                                           wingcol=c("black","black", "gray60" ),
                                           scale.factor = 2,
                                           scale.display = FALSE,
                                           wingframe = FALSE,
                                           winglwd=c(3, 3, 3))
dev.off()

# Ethiopia vs Zambia at 24 deg Males
Eth_Zam_24_Females_dif <-  procoords[7, ] -  procoords[5, ]
Eth_Zam_24_Females_mean <- colMeans(procoords[c(7,5),])
EthZam_24_Females_PD <- norm_vec(Eth_Zam_24_Females_dif)

pdf(file="../outputs/WingEffectPlots/wingeffect_Eth_Zam_t24_Females.pdf")
wingEffect_EthZam_t24_Females <- WingEffect( Eth_Zam_24_Females_mean, Eth_Zam_24_Females_dif , 
                                           Eth_Zam_24_Females_dif,
                                           wingcol=c("black","black", "gray60" ),
                                           scale.factor = 2,
                                           scale.display = FALSE,
                                           wingframe = FALSE,
                                           winglwd=c(3, 3, 3))
dev.off()

# Ethiopia vs Zambia at 28 deg Males
Eth_Zam_28_Females_dif <-  procoords[11, ] -  procoords[9, ]
Eth_Zam_28_Females_mean <- colMeans(procoords[c(11,9),])
EthZam_28_Females_PD <- norm_vec(Eth_Zam_28_Females_dif)

pdf(file="../outputs/WingEffectPlots/wingeffect_Eth_Zam_t28_Females.pdf")
wingEffect_EthZam_t24_Females <- WingEffect( Eth_Zam_28_Females_mean, Eth_Zam_28_Females_dif , 
                                           Eth_Zam_28_Females_dif,
                                           wingcol=c("black","black", "gray60" ),
                                           scale.factor = 2,
                                           scale.display = FALSE,
                                           wingframe = FALSE,
                                           winglwd=c(3, 3, 3))
dev.off()



# making wingplots into ggplot friendly objects
as <- image_read_pdf("../outputs/WingEffectPlots/cropped/wingeffect_Eth_Zam_t18_Females.pdf")

at <- image_read_pdf("../outputs/WingEffectPlots/cropped/wingeffect_Eth_Zam_t24_Females.pdf")

au <- image_read_pdf("../outputs/WingEffectPlots/cropped/wingeffect_Eth_Zam_t28_Females.pdf")

av <- image_read_pdf("../outputs/WingEffectPlots/cropped/wingeffect_Eth_Zam_t18_males.pdf")

aw <- image_read_pdf("../outputs/WingEffectPlots/cropped/wingeffect_Eth_Zam_t24_males.pdf")

ax <- image_read_pdf("../outputs/WingEffectPlots/cropped/wingeffect_Eth_Zam_t28_males.pdf")


# turing wingplot pdf into ggplot friendly object
ay <- ggdraw() + draw_image(as, scale = 1) + 
  annotate("text", x = 0.5, y = c(0.9, 0.08), label = c("18 deg", "PD = 0.011"), size = 24)


az <- ggdraw() + draw_image(at, scale = 1) + 
  annotate("text", x = 0.5, y = c(0.92, 0.08), label = c("24 deg", "PD = 0.014"), size = 24)

ba <- ggdraw() + draw_image(au, scale = 1) + 
  annotate("text", x = 0.5, y = c(0.92, 0.08), label = c("28 deg", "PD = 0.016"), size = 24)


bb <- ggdraw() + draw_image(av, scale = 1) + 
  annotate("text", x = 0.5, y = c(0.92, 0.08), label = c("18 deg",  "PD = 0.012"), size = 24)

bc <- ggdraw() + draw_image(aw, scale = 1) + 
  annotate("text", x = 0.5, y = c(0.92, 0.08), label = c("24 deg", "PD = 0.014"), size = 24)


bd <- ggdraw() + draw_image(ax, scale = 1) + 
  annotate("text", x = 0.5, y = c(0.92, 0.08), label = c("28 deg", "PD = 0.015"), size = 24)

# Putting together the Procurstes Distance between mean Ethiopian and Zambian wing shape at each temperature
EthZam_PD <- rbind(EthZam_18_Males_PD, EthZam_24_Males_PD, EthZam_28_Males_PD, EthZam_18_Females_PD, EthZam_24_Females_PD, EthZam_28_Females_PD)

EthZam_PD_data <- cbind(temp_wings_shape_means[c(2,6,10,1,5,9),c(1,3)], EthZam_PD)

# Calculating the correlations between the mean shape vectors for the ethiopian and zambian populations at different temperatures to look at whether the changes seen are in direction
Ethiopia_18_Male_Vec <- as.numeric(temp_wings_shape_means[4,4:99])
Zambia_18_Male_Vec <- as.numeric(temp_wings_shape_means[2,4:99])

Ethiopia_18_Female_Vec <- as.numeric(temp_wings_shape_means[3,4:99])
Zambia_18_Female_Vec <- as.numeric(temp_wings_shape_means[1,4:99])

Ethiopia_24_Male_Vec <- as.numeric(temp_wings_shape_means[8,4:99])
Zambia_24_Male_Vec <- as.numeric(temp_wings_shape_means[6,4:99])

Ethiopia_24_Female_Vec <- as.numeric(temp_wings_shape_means[7,4:99])
Zambia_24_Female_Vec <- as.numeric(temp_wings_shape_means[5,4:99])

Ethiopia_28_Male_Vec <- as.numeric(temp_wings_shape_means[12,4:99])
Zambia_28_Male_Vec <- as.numeric(temp_wings_shape_means[10,4:99])

Ethiopia_28_Female_Vec <- as.numeric(temp_wings_shape_means[11,4:99])
Zambia_28_Female_Vec <- as.numeric(temp_wings_shape_means[9,4:99])

EthZam_18_Males_cor <- cor(Ethiopia_18_Male_Vec, Zambia_18_Male_Vec)
EthZam_24_Males_cor <- cor(Ethiopia_24_Male_Vec, Zambia_24_Male_Vec)
EthZam_28_Males_cor <- cor(Ethiopia_28_Male_Vec, Zambia_28_Male_Vec)
EthZam_18_Females_cor <- cor(Ethiopia_18_Female_Vec, Zambia_18_Female_Vec)
EthZam_24_Females_cor <- cor(Ethiopia_24_Female_Vec, Zambia_24_Female_Vec)
EthZam_28_Females_cor <- cor(Ethiopia_28_Female_Vec, Zambia_28_Female_Vec)

EthZam_Cor <- rbind(EthZam_18_Males_cor, EthZam_24_Males_cor, EthZam_28_Males_cor, EthZam_18_Females_cor, EthZam_24_Females_cor, EthZam_28_Females_cor)

EthZam_cor_data <- cbind(temp_wings_shape_means[c(2,6,10,1,5,9),c(1,3)], EthZam_Cor)
###############################################################

################# Figure 5 ##################

legend_Figure5<- get_legend(aj_leg + theme(legend.position = "bottom", legend.title = element_text(size = 65), legend.text = element_text(size = 65), legend.key.size = unit(3,"cm")))

Figure5 <- plot_grid(af, ah_sex_corrected, aj, am, axis = "tblr", labels = "AUTO", label_size = 65)

Figure5with_Legend <- plot_grid( Figure5, legend_Figure5, ncol = 1, rel_heights = c(1, 0.1), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/Figure5_ReactionNorms_Wingsize_CV_TotVar_rSDE2.pdf", height = 39,  width = 41)
Figure5with_Legend
dev.off()

#############################################

############### Sup Figure 9 ################

legend_SupFig9<- get_legend(ag_leg + theme(legend.position = "bottom", legend.box = "vertical", 
                                         legend.text = element_text(size = 65), 
                                         legend.key.size = unit(3,"cm")))

SupFig9<- plot_grid(ag, ai_sex_corrected, al, ak, axis = "tblr", labels = "AUTO", label_size = 65)

SupFig9_withlegend <- plot_grid(SupFig9, legend_SupFig9, ncol = 1, rel_heights = c(1, .1), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/SupFigure9_MacroCan_Mut_LevDev_rSDE_Ecc.pdf", height = 39,  width = 41)
SupFig9_withlegend
dev.off()

#############################################


############### SupFigure 10 #################
legend_SupFig10<- get_legend(an_leg + theme(legend.position = "bottom", legend.text = element_text(size = 65),
                                          legend.title = element_text(size = 65),
                                          legend.key.size = unit(3,"cm")))

SupFig10 <- plot_grid(an, ap, aq, ar, axis = "tblr", labels = "AUTO", label_size = 65)

SupFig10_withlegend <- plot_grid(SupFig10, legend_SupFig10, ncol = 1, rel_heights = c(1, .1), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/SupFigure10_MacroCan_CV_totvar_vs_prop_mut.pdf", height = 35,  width = 38)
SupFig10_withlegend
dev.off()
#############################################

########### Supplementary Figure 11: Wing Plots ###############

SupFig11_males <- plot_grid(bb, bc, bd, nrow = 3, axis = "tblr", labels = "AUTO", label_size = 65)

pdf(file = "../outputs/PaperFigures/SupFigure11_TemperatureWingPlots_males.pdf", height = 38,  width = 20)
SupFig11_males
dev.off()


SupFig11_females <- plot_grid(ay, az, ba, nrow = 3, axis = "tblr", labels = "AUTO", label_size = 65)

pdf(file = "../outputs/PaperFigures/SupFigure11_TemperatureWingPlots_females.pdf", height = 38,  width = 20)
SupFig11_females
dev.off()