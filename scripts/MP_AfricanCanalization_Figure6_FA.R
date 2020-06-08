######################################################################
########## Fluctuating Assymertry Figure 6, Supplementary figures 12, 13, and 14 script ###########
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
library(data.table)

###############################################################


######################## Load Data ############################
#FA Data
source('./MP_African_Canalization_FA_cleanup_outlier.R')

###############################################################

#################### Set ggplot theme #########################

pref_theme <- theme_classic() + 
  theme(text = element_text(size = 65), 
        axis.title.y = element_text(margin = margin(t = 10, r = 20, b = 10, l = 10)),
        axis.title.x = element_text(margin = margin(t = 20, r = 10, b = 10, l = 10)),
        axis.text = element_text(margin = margin(t = 15, r = 15, b = 15, l = 15)))
theme_set(pref_theme)

###############################################################

######################### Analysis ############################
wings_SizeFA_NoOut_Temp$population <- relevel(wings_SizeFA_NoOut_Temp$population, "Ethiopia")

#Model and Plot for FA1

wings_SizeFA_NoOut_Temp$abs_diff_transf <- wings_SizeFA_NoOut_Temp$abs_diff + 0.0001


lmm.FA1_gamma <- glmer(abs_diff_transf ~ population*sex*temperature + (1 + temperature + sex|line:population),
                              data = wings_SizeFA_NoOut_Temp, family = "Gamma",
                              control=glmerControl(optimizer="bobyqa",
                                                   optCtrl=list(maxfun=2e5)))


summary(lmm.FA1_gamma)

#Supplementary Table S21
car::Anova(lmm.FA1_gamma)
lmm.FA1_gamma_marginal_means <- as.data.frame(Effect(c("population","sex",  "temperature"), lmm.FA1_gamma))


FA1_lines <- aggregate(wings_SizeFA_NoOut_Temp[,c(12,13,14,15,16,17,18,19)],
                             c(wings_SizeFA_NoOut_Temp["sex"],
                               wings_SizeFA_NoOut_Temp["nutrition"],
                               wings_SizeFA_NoOut_Temp["temperature"],
                               wings_SizeFA_NoOut_Temp["population"],
                               wings_SizeFA_NoOut_Temp["temp_zeroed"],
                               wings_SizeFA_NoOut_Temp["line"]),
                             mean)

FA1_lines$population <- relevel(FA1_lines$population, "Ethiopia")

# Supplementary Figure S12
be <- ggplot(lmm.FA1_gamma_marginal_means, 
             aes(x=temperature, y=fit, colour = population, shape = sex ))+
  geom_point(size = 14, position = position_dodge(0.7), show.legend = FALSE) +
  geom_errorbar(aes(ymin=lower, ymax=upper), lwd = 1.2, width = 0.2, 
                position = position_dodge(0.7), show.legend = FALSE)+
  geom_jitter(data = wings_SizeFA_NoOut_Temp, aes(x = temperature, y = abs_diff_transf), size = 6, 
              position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.9), legend.box = "horizontal")+
  scale_y_continuous(name="Wing Size FA (FA1)")+
  scale_x_discrete(name=expression("Temperature " ( degree*C)),
                   labels = c(18, 24, 28))+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population",
                      labels = c("High-Altitude",
                                 "Low-Altitude"))+
  scale_shape_discrete(name  ="Sex", labels = c("Female","Male"))

be

be_leg <- ggplot(lmm.FA1_gamma_marginal_means, 
             aes(x=temperature, y=fit, colour = population, shape = sex))+
  geom_point(size = 14, position = position_dodge(0.7), show.legend = TRUE) +
  geom_errorbar(aes(ymin=lower, ymax=upper), lwd = 1.2, width = 0.2, 
                position = position_dodge(0.7), show.legend = FALSE)+
  geom_jitter(data = wings_SizeFA_NoOut_Temp, aes(x = temperature, y = abs_diff_transf), size = 6, 
              position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = "bottom", legend.box = "vertical")+
  scale_y_continuous(name="Wing Size FA (FA1)")+
  scale_x_discrete(name=expression("Temperature " ( degree*C)),
                   labels = c(18, 24, 28))+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population",
                      labels = c("High-Altitude",
                                 "Low-Altitude"))+
  scale_shape_discrete(name  ="Sex", labels = c("Female","Male"))

be_leg


# Model and plot for FA8

wings_SizeFA_NoOut_Temp$abs_ln_ratio_transf <- wings_SizeFA_NoOut_Temp$abs_ln_ratio + 0.0001


FA8_lines <- aggregate(wings_SizeFA_NoOut_Temp["abs_ln_ratio_transf"],
                       c(wings_SizeFA_NoOut_Temp["temperature"],
                         wings_SizeFA_NoOut_Temp["population"],
                         wings_SizeFA_NoOut_Temp["sex"],
                         wings_SizeFA_NoOut_Temp["temp_zeroed"],
                         wings_SizeFA_NoOut_Temp["line"]),
                       mean, na.rm = FALSE)

FA8_lines$population <- relevel(FA8_lines$population, "Ethiopia")


lmm.FA8_gamma <- glmer(abs_ln_ratio_transf ~ sex*population*temperature + 
                         (1 + temperature + sex|line:population),
                       data = wings_SizeFA_NoOut_Temp, family = "Gamma",
                       control=glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))


summary(lmm.FA8_gamma)

# Supplementary Table S22
car::Anova(lmm.FA8_gamma)
lmm.FA8_gamma_marginal_means <- as.data.frame(Effect(c("sex", "population", "temperature"), lmm.FA8_gamma))


# Figure 6A
bf <- ggplot(lmm.FA8_gamma_marginal_means, 
             aes(x=temperature, y=fit, colour = population, shape = sex))+
  geom_point(size = 14, position = position_dodge(0.7), show.legend = FALSE) +
  geom_errorbar(aes(ymin=lower, ymax=upper), lwd = 1.2, width = 0.2, 
                position = position_dodge(0.7), show.legend = FALSE)+
  geom_jitter(data = wings_SizeFA_NoOut_Temp, aes(x = temperature, y = abs_ln_ratio_transf), size = 6, 
              position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.9), legend.box = "horizontal")+
  scale_y_continuous(name="Wing Size FA (FA8)")+
  scale_x_discrete(name=expression("Temperature " ( degree*C)),
                   labels = c(18, 24, 28))+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))+
  scale_shape_discrete(name  ="Sex", labels = c("Female","Male"))

bf

# Model and Plot for PD between left and right

PD_left_right$population <- relevel(PD_left_right$population, "Ethiopia")

lmm1.PD_LR_gamma <- glmer(PD_left_right ~ sex*population*temperature + 
                         (1 + sex + temperature|line:population),
                       data = subset(PD_left_right, nutrition == "f100"), family = "Gamma",
                       control=glmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))

summary(lmm1.PD_LR_gamma)

# Supplementary Table S23
car::Anova(lmm1.PD_LR_gamma)

lmm1.PD_LR_gamma_marginal_means <- as.data.frame(Effect(c("sex", "population", "temperature"), lmm1.PD_LR_gamma))


PD_left_right_temp <- subset(PD_left_right, nutrition == "f100")


PD_left_right_temp_lines <- aggregate(PD_left_right_temp["PD_left_right"],
                                      by = c(PD_left_right_temp["line"],
                                             PD_left_right_temp["sex"],
                                             PD_left_right_temp["temperature"],
                                             PD_left_right_temp["population"]),
                                      mean,  na.rm = TRUE)


PD_left_right_temp_lines$population <- relevel(PD_left_right_temp_lines$population, "Ethiopia")

# Figure S11B
bg<- ggplot(lmm1.PD_LR_gamma_marginal_means, 
            aes(x=temperature, y=fit, colour = population, shape = sex))+
  geom_point(size = 14, position = position_dodge(0.7), show.legend = FALSE) +
  geom_errorbar(aes(ymin=lower, ymax=upper), lwd = 1.5, width = 0.3, 
                position = position_dodge(0.7), show.legend = FALSE)+
  geom_jitter(data = PD_left_right_temp, aes(x = temperature, y = PD_left_right), size = 6, position = position_jitterdodge(jitter.width = 0.2), show.legend = FALSE, alpha = 0.5)+
  theme(legend.position = c(0.5, 0.9), legend.box = "horizontal")+
  scale_y_continuous(name = expression("Wing Shape FA" (PD["LR"])))+
  scale_x_discrete(name = expression("Temperature " ( degree*C)),
                   labels = c(18, 24, 28))+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))+
  scale_shape_discrete(name  ="Sex", labels = c("Female","Male"))

bg


############################################################

#################### Measurement Error #####################

# Wing size ME
RepMes_only <- merge(wings_size_shape_FA, wings_ME, by = c("line",
                                                           "condition",
                                                           "replicate",
                                                           "sex",
                                                           "individual",
                                                           "nutrition",
                                                           "temperature",
                                                           "temp_num",
                                                           "temp_zeroed",
                                                           "population",
                                                           "side",
                                                           "inds"))


RepMes_only_size <- RepMes_only[,c(1:13,112,14,213)]
RepMes_only_size_1 <- RepMes_only_size[-c(43,94,109,112),]

RepMes_only_size_long <- gather(RepMes_only_size_1, trait, measurement, CS.x: CS.y, factor_key=TRUE)


lm_ME <- lm(measurement ~ side + inds + side:inds , data = RepMes_only_size_long)

#Supplementary Table S24

anova(lm_ME)
#stargazer::stargazer(anova(lm_ME), summary = FALSE)

summary(lm_ME)


# Wing Shape ME
RepMes_shapesize_Mes1 <- RepMes_only[,c(1:12, 13, 14, 15:110)]
RepMes_shapesize_Mes2 <- RepMes_only[,c(1:12, 112, 213, 117:212)]

names(RepMes_shapesize_Mes1) <- c("line", "condition", "replicate", "sex", "individual", "nutrition", "temperature", "temp_num",   "temp_zeroed", "population",  "side", "inds", "mesID", "CS", "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "x5", "y5", "x6", "y6", "x7", "y7", "x8", "y8", "x9", "y9", "x10", "y10", "x11",  "y11", "x12", "y12", "x13", "y13", "x14", "y14", "x15", "y15", "x16",  "y16", "x17", "y17", "x18", "y18", "x19", "y19", "x20", "y20", "x21", "y21", "x22", "y22", "x23", "y23", "x24", "y24", "x25", "y25", "x26", "y26", "x27", "y27", "x28", "y28", "x29", "y29", "x30", "y30", "x31", "y31", "x32", "y32", "x33", "y33", "x34", "y34", "x35", "y35", "x36", "y36", "x37", "y37", "x38", "y38", "x39", "y39", "x40",  "y40", "x41", "y41", "x42", "y42", "x43", "y43", "x44", "y44", "x45", "y45", "x46", "y46", "x47", "y47", "x48", "y48")


names(RepMes_shapesize_Mes2) <- c("line", "condition", "replicate", "sex", "individual", "nutrition", "temperature", "temp_num",   "temp_zeroed", "population",  "side", "inds", "mesID", "CS", "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "x5", "y5", "x6", "y6", "x7", "y7", "x8", "y8", "x9", "y9", "x10", "y10", "x11",  "y11", "x12", "y12", "x13", "y13", "x14", "y14", "x15", "y15", "x16",  "y16", "x17", "y17", "x18", "y18", "x19", "y19", "x20", "y20", "x21", "y21", "x22", "y22", "x23", "y23", "x24", "y24", "x25", "y25", "x26", "y26", "x27", "y27", "x28", "y28", "x29", "y29", "x30", "y30", "x31", "y31", "x32", "y32", "x33", "y33", "x34", "y34", "x35", "y35", "x36", "y36", "x37", "y37", "x38", "y38", "x39", "y39", "x40",  "y40", "x41", "y41", "x42", "y42", "x43", "y43", "x44", "y44", "x45", "y45", "x46", "y46", "x47", "y47", "x48", "y48")


RepMes_Shape <- rbind(RepMes_shapesize_Mes1, RepMes_shapesize_Mes2)
RepMes_Shape_1 <- RepMes_Shape[-c(43, 94, 109, 112, 201, 252, 267, 270),]

RepMes_Shape_1$inds <- droplevels(RepMes_Shape_1$inds)

RepMes_Shape_1$indsnum <- as.numeric(as.factor(RepMes_Shape_1$inds))

RepMes_Shape_1$rep <- as.numeric(as.factor(RepMes_Shape_1$mesID))

RepMes_Shape_1$sidenum <- as.numeric(as.factor(RepMes_Shape_1$side))


wings_land<- RepMes_Shape_1[,15:110]
wings_descr <- RepMes_Shape_1[,c(1:13,111, 112, 113)]
wings_csize <- RepMes_Shape_1$CS

#converting landmarks into 3D array because geom
wings_land_3D <- arrayspecs(wings_land, 48, 2)

#converting to geomorph data frame
gdf <- geomorph.data.frame(shape = wings_land_3D,
                           size = wings_csize,
                           side = wings_descr$side,
                           ind = wings_descr$inds,
                           replicate = wings_descr$rep)


WingsME_Sym <- bilat.symmetry(shape, ind = ind , side = side, replicate = replicate, object.sym = FALSE, iter = 100, data = gdf, RRPP = TRUE)

summary(WingsME_Sym)
WingsME_Sym$shape.anova
bh <- plot(WingsME_Sym, warpgrids = TRUE)

############################################################

########################## using bilat.symmetry #################

# for all groups together
FA_shape <- arrayspecs(wings_size_shape_FA[,14:109], 48, 2)
FA_size <- wings_size_shape_FA$CS
FA_side <- as.numeric(as.factor(wings_size_shape_FA$side))
FA_ind <- as.numeric(as.factor(wings_size_shape_FA$inds))
FA_line <- as.factor(wings_size_shape_FA$line)
FA_nutrition <- as.factor(wings_size_shape_FA$nutrition)
FA_temperature <- as.factor(wings_size_shape_FA$temperature)
FA_temp_nut <- as.factor(wings_size_shape_FA$temp_num)
FA_population <- as.factor(wings_size_shape_FA$population)
FA_sex <- as.factor(wings_size_shape_FA$sex)

gdf <- geomorph.data.frame(shape = FA_shape,
                           size = FA_size,
                           side = FA_side,
                           ind = FA_ind,
                           line = FA_line,
                           nutrition = FA_nutrition,
                           temperature = FA_temperature,
                           temp_num = FA_temp_nut,
                           population = FA_population,
                           sex = FA_sex)

# Assuming all groups (sexes, pops ...) have the same DA and AA.
Wings_Sym <- bilat.symmetry(shape, ind = ind , side = side, object.sym = FALSE, iter = 100, data = gdf, RRPP = TRUE)

summary(Wings_Sym)
Wings_Sym$shape.anova

bh <- plot(Wings_Sym, warpgrids = TRUE)

Wings_Sym$asymm.shape
Wings_Sym$symm.shape


DA_shapecomponent <- Wings_Sym$DA.component

FA_shapecomponent <- Wings_Sym$FA.component

dummy_num <- rep(1, 728)

FA_labels <- unique(wings_size_shape_FA$inds)

FA_labels_df <- data.frame(dummy_num, FA_labels)

FA_labels_df$ind <- FA_labels_df$FA_labels

FA_labels_df_1 <- separate(data = FA_labels_df, col = FA_labels, into= c("line", "condition", "replicate", "sex", "individual"), sep = "\\.")

FA_labels_df_1[,2:6] <- lapply(FA_labels_df_1[,2:6], factor)

nutrition <- FA_labels_df_1$condition

levels(nutrition) <- c("f15","f100", "f100", "f100")

temperature <- FA_labels_df_1$condition

levels(temperature) <- c("t24","t18", "t24", "t28")

FA_labels_df_2 <- cbind(FA_labels_df_1, nutrition, temperature)

pop_dummy <- substr(FA_labels_df_2$line, 1, 1)

FA_labels_df_2$population <- ifelse(pop_dummy == "z", "Zambia", "Ethiopia")

FA_labels_df_2$population <- factor(FA_labels_df_2$population)


FA_df_CS <- merge(FA_labels_df_2, wings_SizeFA_3.5, by = c("line", "condition", "replicate", "sex", "individual", "nutrition", "temperature", "population"))


gdf_fa <- geomorph.data.frame(shape = FA_shapecomponent,
                           line = FA_df_CS$line,
                           nutrition = FA_df_CS$nutrition,
                           temperature = FA_df_CS$temperature,
                           population = FA_df_CS$population,
                           sex = FA_df_CS$sex,
                           size = FA_df_CS$meanCS)


# Model testing the effects of wing size, sex, population and temperarure on shape FA after removing DA
FA_model_2 <- procD.lm(shape ~ (size + sex + population   + temperature)^3 + population/line, 
                       data = gdf_fa, iter = 100, subset = nutrition == "f100",
                       RRPP=TRUE, print.progress = TRUE)

summary(FA_model_2)


#Morphological disparity calculating the PD between different groups
FA_model_2_disp <- morphol.disparity(shape ~ size*sex , groups = ~ population*temperature,  partial = FALSE,
                                     iter = 1000, data = gdf_fa, print.progress = TRUE,
                                     subset = nutrition == "f100", RRPP=TRUE)


summary(FA_model_2_disp)


### Calculating the L2 norm of FA component for each individual
FA_shapecomponent_2D <- two.d.array(FA_shapecomponent)

# norm FA component for Zambia
norm_vec <- function(x) sqrt(sum(x^2))

Ind_FA_norm <- apply(FA_shapecomponent_2D, 1, norm_vec)

FA_df_CS$Ind_FA_norm <- Ind_FA_norm

FA_df_CS$Ind_FA_norm_alt <- (FA_df_CS$Ind_FA_norm - 1)


lmm1.FA_PD <- glm(Ind_FA_norm  ~ (population + sex + temperature + meanCS)^3, family = "Gamma",
                   data = subset(FA_df_CS, nutrition == "f100"))


summary(lmm1.FA_PD)
Anova(lmm1.FA_PD)

check_FA_size <- ggplot(data = FA_df_CS, aes(y = Ind_FA_norm_alt, x = meanCS))+
  geom_point()
check_FA_size

cor.test(FA_df_CS$Ind_FA_norm_alt, FA_df_CS$meanCS)


lmm1.FA_PD_nosize_gamma_marginal_means <- as.data.frame(Effect(c("population","sex", "temperature"), lmm1.FA_PD))

zz <- ggplot(lmm1.FA_PD_nosize_gamma_marginal_means, 
             aes(x=temperature, y=(fit-1)*1000, 
                 colour = population, shape = sex ))+
  geom_point(size = 14, position = position_dodge(0.7), 
             show.legend = FALSE) +
  geom_errorbar(aes(ymin=(lower-1)*1000, ymax=(upper-1)*1000 ), 
                lwd = 1.5, width = 0.3, 
                position = position_dodge(0.7), show.legend = FALSE)+
  geom_jitter(data = subset(FA_df_CS, nutrition == "f100"), 
              aes(x = temperature, y = ((Ind_FA_norm-1)*1000) ),
              size = 6, position = position_jitterdodge(jitter.width = 0.2), 
              show.legend = FALSE, alpha = 0.3)+
  theme(legend.position = c(0.5, 0.9), legend.box = "horizontal")+
  scale_y_continuous(name="Wing Shape FA")+
  scale_x_discrete(name=expression("Temperature " ( degree*C)),
                   labels = c(18, 24, 28))+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))+
  scale_shape_discrete(name  ="Sex", labels = c("Female","Male"))

zz


##############################################
### Redoing the analysis using bilat.symmetry separately for each population in order to avoid assuming similar DA
############################################

# Ethiopia
wings_size_shape_FA_Ethiopia <- subset(wings_size_shape_FA, population == "Ethiopia")


FA_shape_Ethiopia <- arrayspecs(wings_size_shape_FA_Ethiopia[,14:109], 48, 2)
FA_size_Ethiopia <- wings_size_shape_FA_Ethiopia$CS
FA_side_Ethiopia <- wings_size_shape_FA_Ethiopia$side
FA_ind_Ethiopia <- wings_size_shape_FA_Ethiopia$inds


gdf_Ethiopia <- geomorph.data.frame(shape = FA_shape_Ethiopia,
                           size = FA_size_Ethiopia,
                           side = FA_side_Ethiopia,
                           ind = FA_ind_Ethiopia)


Wings_Sym_Ethiopia <- bilat.symmetry(shape, ind = ind , side = side,object.sym = FALSE, iter = 1, data = gdf_Ethiopia, RRPP = TRUE)

summary(Wings_Sym_Ethiopia)
Wings_Sym_Ethiopia$shape.anova

plot(Wings_Sym_Ethiopia, warpgrids = TRUE)

FA_shapecomponent_Ethiopia <- Wings_Sym_Ethiopia$FA.component

dummy_num <- rep(1, 302)

wings_SizeFA_3.5_Ethiopia <- subset(wings_SizeFA_3.5, population == "Ethiopia")


FA_labels_Ethiopia <- unique(wings_size_shape_FA_Ethiopia$inds)

FA_labels_df_Ethiopia <- data.frame(dummy_num, FA_labels_Ethiopia)

FA_labels_df_Ethiopia$ind <- FA_labels_df_Ethiopia$FA_labels_Ethiopia

FA_labels_df_1_Ethiopia <- separate(data = FA_labels_df_Ethiopia, col = FA_labels_Ethiopia, into= c("line", "condition", "replicate", "sex", "individual"), sep = "\\.")

FA_labels_df_1_Ethiopia[,2:6] <- lapply(FA_labels_df_1_Ethiopia[,2:6], factor)

nutrition <- FA_labels_df_1_Ethiopia$condition

levels(nutrition) <- c("f15","f100", "f100", "f100")

temperature <- FA_labels_df_1_Ethiopia$condition

levels(temperature) <- c("t24","t18", "t24", "t28")

FA_labels_df_2_Ethiopia <- cbind(FA_labels_df_1_Ethiopia, nutrition, temperature)

pop_dummy <- substr(FA_labels_df_2_Ethiopia$line, 1, 1)

FA_labels_df_2_Ethiopia$population <- ifelse(pop_dummy == "z", "Zambia", "Ethiopia")

FA_labels_df_2_Ethiopia$population <- factor(FA_labels_df_2_Ethiopia$population)

FA_df_CS_Ethiopia <- merge(FA_labels_df_2_Ethiopia, wings_SizeFA_3.5_Ethiopia, by = c("line", "condition", "replicate", "sex", "individual", "nutrition", "temperature", "population"))

gdf_Ethiopia_fa <- geomorph.data.frame(shape = FA_shapecomponent_Ethiopia,
                              size = FA_df_CS_Ethiopia$meanCS,
                              line = FA_df_CS_Ethiopia$line,
                              nutrition = FA_df_CS_Ethiopia$nutrition,
                              temperature = FA_df_CS_Ethiopia$temperature,
                              sex = FA_df_CS_Ethiopia$sex,
                              population = FA_df_CS_Ethiopia$population,
                              ind = FA_df_CS_Ethiopia$ind)


FA_model_1_Ethiopia <- procD.lm(shape ~ size*sex*temperature, 
                       data = gdf_Ethiopia_fa, iter = 100, subset = nutrition == "f100",
                       RRPP=TRUE, print.progress = TRUE)


summary(FA_model_1_Ethiopia)


FA_model_1_disp_Ethiopia <- morphol.disparity(shape ~ size , groups = ~ population*sex*temperature,  
                                     iter = 5000, data = gdf_Ethiopia_fa, print.progress = TRUE, 
                                     subset = nutrition == "f100", RRPP=TRUE)



### Calculating the L2 norm of FA component for each individual
FA_shapecomponent_2D_Ethiopia <- two.d.array(FA_shapecomponent_Ethiopia)

# norm FA component for Zambia
norm_vec <- function(x) sqrt(sum(x^2))

Ind_FA_norm_Ethiopia <- apply(FA_shapecomponent_2D_Ethiopia, 1, norm_vec)

FA_df_CS_Ethiopia$Ind_FA_norm <- Ind_FA_norm_Ethiopia


lmm1.FA_PD_Ethiopia <- lmer(Ind_FA_norm ~ sex*temperature + 
                     (1 |line:population),
                   data = subset(FA_df_CS_Ethiopia, nutrition == "f100"))

summary(lmm1.FA_PD_Ethiopia)


car::Anova(lmm1.FA_PD_Ethiopia)

lmm1.FA_PD_Ethiopia_marginal_means <- as.data.frame(Effect(c("sex", "temperature"), lmm1.FA_PD_Ethiopia))


# Zambia

wings_size_shape_FA_Zambia <- subset(wings_size_shape_FA, population == "Zambia")


FA_shape_Zambia <- arrayspecs(wings_size_shape_FA_Zambia[,14:109], 48, 2)
FA_size_Zambia <- wings_size_shape_FA_Zambia$CS
FA_side_Zambia <- as.numeric(as.factor(wings_size_shape_FA_Zambia$side))
FA_ind_Zambia <- as.numeric(as.factor(wings_size_shape_FA_Zambia$inds))


gdf_Zambia <- geomorph.data.frame(shape = FA_shape_Zambia,
                                    size = FA_size_Zambia,
                                    side = FA_side_Zambia,
                                    ind = FA_ind_Zambia)


Wings_Sym_Zambia <- bilat.symmetry(shape, ind = ind , side = side,object.sym = FALSE, iter = 1, data = gdf_Zambia, RRPP = TRUE)

summary(Wings_Sym_Zambia)
Wings_Sym_Zambia$shape.anova

plot(Wings_Sym_Zambia, warpgrids = TRUE)


FA_shapecomponent_Zambia <- Wings_Sym_Zambia$FA.component


dummy_num_zam <- rep(1, 426)

wings_SizeFA_3.5_Zambia <- subset(wings_SizeFA_3.5, population == "Zambia")


FA_labels_Zambia <- unique(wings_size_shape_FA_Zambia$inds)

FA_labels_df_Zambia <- data.frame(dummy_num_zam, FA_labels_Zambia)

FA_labels_df_Zambia$ind <- FA_labels_df_Zambia$FA_labels_Zambia


FA_labels_df_1_Zambia <- separate(data = FA_labels_df_Zambia, col = FA_labels_Zambia, into= c("line", "condition", "replicate", "sex", "individual"), sep = "\\.")

FA_labels_df_1_Zambia[,2:6] <- lapply(FA_labels_df_1_Zambia[,2:6], factor)

nutrition <- FA_labels_df_1_Zambia$condition

levels(nutrition) <- c("f15","f100", "f100", "f100")

temperature <- FA_labels_df_1_Zambia$condition

levels(temperature) <- c("t24","t18", "t24", "t28")

FA_labels_df_2_Zambia <- cbind(FA_labels_df_1_Zambia, nutrition, temperature)

pop_dummy <- substr(FA_labels_df_2_Zambia$line, 1, 1)

FA_labels_df_2_Zambia$population <- ifelse(pop_dummy == "z", "Zambia", "Ethiopia")

FA_labels_df_2_Zambia$population <- factor(FA_labels_df_2_Zambia$population)

FA_df_CS_Zambia <- merge(FA_labels_df_2_Zambia, wings_SizeFA_3.5_Zambia, by = c("line", "condition", "replicate", "sex", "individual", "nutrition", "temperature", "population"))

gdf_Zambia_fa <- geomorph.data.frame(shape = FA_shapecomponent_Zambia,
                                       size = FA_df_CS_Zambia$meanCS,
                                       line = FA_df_CS_Zambia$line,
                                       nutrition = FA_df_CS_Zambia$nutrition,
                                       temperature = FA_df_CS_Zambia$temperature,
                                       sex = FA_df_CS_Zambia$sex,
                                       population = FA_df_CS_Zambia$population)


FA_model_1_Zambia <- procD.lm(shape ~ size*sex*temperature, 
                                data = gdf_Zambia_fa, iter = 100, subset = nutrition == "f100",
                                RRPP=TRUE, print.progress = TRUE)


summary(FA_model_1_Zambia)


FA_model_1_disp_Zambia <- morphol.disparity(shape ~ size , groups = ~ population*sex*temperature,  
                                              iter = 5000, data = gdf_Zambia_fa, print.progress = TRUE, 
                                              subset = temperature == "t24", RRPP=TRUE)


### Calculating the L2 norm of FA component for each individual
FA_shapecomponent_2D_Zambia <- two.d.array(FA_shapecomponent_Zambia)

# norm FA component for Zambia
norm_vec <- function(x) sqrt(sum(x^2))

Ind_FA_norm_Zambia <- apply(FA_shapecomponent_2D_Zambia, 1, norm_vec)

FA_df_CS_Zambia$Ind_FA_norm <- Ind_FA_norm_Zambia


lmm1.FA_PD_Zambia <- lmer(Ind_FA_norm ~ sex*temperature + 
                              (1 |line:population),
                            data = subset(FA_df_CS_Zambia, nutrition == "f100"))

summary(lmm1.FA_PD_Zambia)


car::Anova(lmm1.FA_PD_Zambia)

lmm1.FA_PD_Zambia_marginal_means <- as.data.frame(Effect(c("sex", "temperature"), lmm1.FA_PD_Zambia))


#################################################################

##Comparing FA with measures of variablity for wing size and shape##

source('MP_AfricanCanalization_Figure5_Macroenvironmental_Canalization.R')

# wing size varaition measures vs FA measures

cv_dat2_fa_only <- cv_dat2[c(6,8,16,17,28,30,37, 38, 48, 49, 56, 57), ]

FA1_lines_comp <- FA1_lines[, c(1,3,4,6,7,8,9,10,11,14)]

FA1_lines_comp_meansex <- aggregate(FA1_lines_comp[,5:10],
                                    c(FA1_lines_comp["temperature"],
                                      FA1_lines_comp["population"],
                                      FA1_lines_comp["line"]),
                                    mean, na.rm = TRUE)


FA8_lines_comp <- FA8_lines[, c(1, 2, 3, 5, 6)]

FA8_lines_comp_meansex <- aggregate(FA8_lines_comp["abs_ln_ratio_transf"],
                                    c(FA1_lines_comp["temperature"],
                                      FA1_lines_comp["population"],
                                      FA1_lines_comp["line"]),
                                    mean, na.rm = TRUE)

PD_left_right_lines_meansex <- aggregate(PD_left_right_temp_lines["PD_left_right"],
                                      c(PD_left_right_temp_lines["temperature"],
                                        PD_left_right_temp_lines["population"],
                                        PD_left_right_temp_lines["line"]),
                                      mean, na.rm = TRUE)

all_FA_size_mes <- merge(FA1_lines_comp_meansex, FA8_lines_comp_meansex)

all_FA_size_shape <- merge(all_FA_size_mes, PD_left_right_lines_meansex)

cv_FA_size_shape_dat <- merge(cv_dat2_fa_only, all_FA_size_shape)

##########################################################################
### Plots for FA analyses vs measures of wing size and shape variation ###
##########################################################################

zza <- ggplot(data = cv_FA_size_shape_dat, aes(x = cv_out, y = abs_diff_transf, 
                                               colour = population, 
                                               shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA1)")+
  scale_x_continuous(name="Wing size CV")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zza

zza_leg <- ggplot(data = cv_FA_size_shape_dat, aes(x = cv_out, y = abs_diff_transf, 
                                               colour = population, 
                                               shape = temperature)) +
  geom_point(size = 14)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA1)")+
  scale_x_continuous(name="Wing size CV")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))+
  scale_shape_discrete(name = "Temperature",
                       labels = c("18 deg", "24 deg", "28 deg"))
zza_leg

zzb <- ggplot(data = cv_FA_size_shape_dat, aes(x = cv_out, y = abs_ln_ratio_transf, 
                                               colour = population, 
                                               shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA8)")+
  scale_x_continuous(name="Wing size CV")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzb


zzc <- ggplot(data = cv_FA_size_shape_dat, aes(x = cv_out, y = PD_left_right, 
                                               colour = population, 
                                               shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing shape FA")+
  scale_x_continuous(name="Wing size CV")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzc


zzd <- ggplot(data = cv_FA_size_shape_dat, aes(x = ls1_out, y = abs_diff_transf, 
                                               colour = population, 
                                               shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA1)")+
  scale_x_continuous(name="Wing size LD")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzd

zze <- ggplot(data = cv_FA_size_shape_dat, aes(x = ls1_out, y = abs_ln_ratio_transf, 
                                               colour = population, 
                                               shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA8)")+
  scale_x_continuous(name="Wing size LD")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zze

zzf <- ggplot(data = cv_FA_size_shape_dat, aes(x = ls1_out, y = PD_left_right, 
                                               colour = population, 
                                               shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing shape FA")+
  scale_x_continuous(name="Wing size LD")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzf

# comparing with wing shape measures of variability

integration_measures_temp_fa_lines_only <- integration_measures_temp[c(6,8,16,17,28,30,37,38,48, 49, 56, 57),]

int_FA_size_shape_dat <- merge(integration_measures_temp_fa_lines_only, all_FA_size_shape)


zzg <- ggplot(data = int_FA_size_shape_dat, aes(x = tot_var_temp_alt, 
                                                y = abs_diff_transf, 
                                               colour = population, 
                                               shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA1)")+
  scale_x_continuous(name="Wing shape TV")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzg

zzh <- ggplot(data = int_FA_size_shape_dat, aes(x = tot_var_temp_alt, 
                                                y = abs_ln_ratio_transf, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA8)")+
  scale_x_continuous(name="Wing shape TV")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzh

zzi <- ggplot(data = int_FA_size_shape_dat, aes(x = tot_var_temp_alt, 
                                                y = PD_left_right, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing shape FA")+
  scale_x_continuous(name="Wing shape TV")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzi


zzj <- ggplot(data = int_FA_size_shape_dat, aes(x = eccentricity_temp, 
                                                y = abs_diff_transf, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA1)")+
  scale_x_continuous(name="Wing shape Ecc")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzj

zzk <- ggplot(data = int_FA_size_shape_dat, aes(x = eccentricity_temp, 
                                                y = abs_ln_ratio_transf, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA8)")+
  scale_x_continuous(name="Wing shape Ecc")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzk

zzl <- ggplot(data = int_FA_size_shape_dat, aes(x = eccentricity_temp, 
                                                y = PD_left_right, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing shape FA")+
  scale_x_continuous(name="Wing shape Ecc")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzl


zzm <- ggplot(data = int_FA_size_shape_dat, aes(x = rSDE_temp_alt, 
                                                y = abs_diff_transf, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA1)")+
  scale_x_continuous(name="Wing shape rSDE")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzm

zzn <- ggplot(data = int_FA_size_shape_dat, aes(x = rSDE_temp_alt, 
                                                y = abs_ln_ratio_transf, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA8)")+
  scale_x_continuous(name="Wing shape rSDE")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzn

zzo <- ggplot(data = int_FA_size_shape_dat, aes(x = rSDE_temp_alt, 
                                                y = PD_left_right, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing shape FA")+
  scale_x_continuous(name="Wing shape rSDE")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzo

zzp <- ggplot(data = int_FA_size_shape_dat, aes(x = rSDE2_temp_alt, 
                                                y = abs_diff_transf, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA1)")+
  scale_x_continuous(name="Wing shape rSDE2")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzp

zzq <- ggplot(data = int_FA_size_shape_dat, aes(x = rSDE2_temp_alt, 
                                                y = abs_ln_ratio_transf, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing size FA (FA8)")+
  scale_x_continuous(name="Wing shape rSDE2")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzq

zzr <- ggplot(data = int_FA_size_shape_dat, aes(x = rSDE2_temp_alt, 
                                                y = PD_left_right, 
                                                colour = population, 
                                                shape = temperature)) +
  geom_point(size = 14, show.legend = FALSE)+
  theme(legend.position = "bottom", legend.box = "horizontal")+
  scale_y_continuous(name="Wing shape FA")+
  scale_x_continuous(name="Wing shape rSDE2")+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude"))
zzr

#################################################################


########### Figure 6 ################


legend_Figure6<- get_legend(be_leg + theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 50), legend.title = element_text(size = 60),  legend.key.size = unit(2,"cm")))

Figure6 <- plot_grid(bf, zz, axis = "tblr", labels = "AUTO", label_size = 65)

Figure6with_Legend <- plot_grid( Figure6, legend_Figure6, ncol = 1, rel_heights = c(1, 0.2), axis = "t", align = "v")


pdf(file = "../outputs/PaperFigures/Figure6_FA8_PD.pdf", height = 18,  width = 32)
Figure6with_Legend
dev.off()

#####################################

########### SupFigure 12 ############


SupFigure12 <- plot_grid(be, bg, axis = "tblr", labels = "AUTO", label_size = 65)
SupFigure12_Legend <- plot_grid( SupFigure12, legend_Figure6, ncol = 1, rel_heights = c(1, 0.2), axis = "t", align = "v")

pdf(file = "../outputs/PaperFigures/SupFigure12_FA1_PD_LR.pdf", height = 18,  width = 32)
SupFigure12_Legend
dev.off()

#########################################

############## SupFigure 13 #############
legend_SupFig13<- get_legend(zza_leg + theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 60), legend.title = element_text(size = 60),  legend.key.size = unit(2,"cm")))

SupFigure13 <- plot_grid(zza, zzb, zzc, zzd, zze, zzf, axis = "tblr", labels = "AUTO", label_size = 65)

SupFigure13_Legend <- plot_grid( SupFigure13, legend_SupFig13, ncol = 1, rel_heights = c(1, 0.1), axis = "t", align = "v")


pdf(file = "../outputs/PaperFigures/SupFigure13_FA_vs_wingsizeVar.pdf", height = 30,  width = 40)
SupFigure13_Legend
dev.off()


############## SupFigure 14 #############

SupFigure14 <- plot_grid(zzg, zzh, zzi, zzj, zzk, zzl, zzm, zzn, zzo, zzp, zzq, zzr, axis = "tblr", nrow = 4, labels = "AUTO", label_size = 65)

SupFigure14_Legend <- plot_grid( SupFigure14, legend_SupFig13, ncol = 1, rel_heights = c(1, 0.1), axis = "t", align = "v")


pdf(file = "../outputs/PaperFigures/SupFigure14_FA_vs_wingshapeVar.pdf", height = 60,  width = 60)
SupFigure14_Legend
dev.off()
