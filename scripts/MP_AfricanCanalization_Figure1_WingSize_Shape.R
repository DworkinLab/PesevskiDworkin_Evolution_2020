###############################################################
###### Wing size and shape analysis and Figure 1 script #######
###############################################################

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
library(pdftools)

source('../misc/WINGPLOTSOURCE.R', chdir = TRUE)

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

################### Wing size raw data ########################

#FIGURE 1A 
a <- ggplot(wings, aes(x=population:sex, y=CS_rescaled, color = population, shape = population)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, lwd=2) + 
  labs(y = "Wing size", x = "") +
  geom_jitter(size = 12, alpha = 0.4, position=position_jitter(0.3), show.legend = TRUE) +
  theme(legend.position = c(0.9,0.9))+
  scale_colour_manual(name = "Population",
                      labels = c("High-Altitude", "Low-Altitude"),
                      values = c("black","gray40")) +
  scale_x_discrete(labels=c("High-Altitude \nFemale","High-Altitude \nMale", 
                            "Low-altitude \nFemale", "Low-Altitude \nMale"))+
  scale_shape_discrete(name = "Population",
                       labels = c("High-Altitude", "Low-Altitude"))

a

###############################################################

###### Microenvironmental Canalization Wing size Analysis ######

#Wing size model
lmm.size <- lmer(CS_rescaled ~ sex*population + (1 + sex|line:population), data = wings)
summary(lmm.size)

#Table S4
Anova(lmm.size)

#Extracting effects from model
lmm.size_marginal_means <- as.data.frame(Effect(c("sex", "population"), lmm.size))

###############################################################


######## Landmarks and semilandmarks ##########
# Figure 1B - Landmarks and semilandmarks
b_1 <- image_read("../outputs/Figure1B_Landmarks.png")

# turning landmark picture into ggplot friendly object
b <- ggdraw() + draw_image(b_1, scale = 1.1)
###############################################################


############## Wing shape raw data - Wing plot ################
#subsetting data for each group
EthFemale_data <- subset(wings, subset = sex == "f" & population == "Ethiopia", drop = T )
EthMale_data <- subset(wings, subset = sex == "m" & population == "Ethiopia", drop = T )
ZamFemale_data <- subset(wings, subset = sex == "f" & population == "Zambia", drop = T )
ZamMale_data <- subset(wings, subset = sex == "m" & population == "Zambia", drop = T )

#making shape matrices 
EthFemale_vec <- data.matrix(EthFemale_data[,7:102])

EthMale_vec <- data.matrix(EthMale_data[,7:102])

ZamFemale_data <- data.matrix(ZamFemale_data[,7:102])

ZamMale_data <- data.matrix(ZamMale_data[,7:102])


#calculating the distance vector between the mean shapes for each pop sex
EthZam_Fem_Diff_Vec <- colMeans(EthFemale_vec) - colMeans(ZamFemale_data)

EthZam_Mal_Diff_Vec <- colMeans(EthMale_vec) - colMeans(ZamMale_data)

#calculating the magnitude of the distance vector aka PD
norm_vec <- function(x) sqrt(sum(x^2))

norm_vec(as.matrix(EthZam_Fem_Diff_Vec))
norm_vec(as.matrix(EthZam_Mal_Diff_Vec))


wings_id <- aggregate( wings[,10:105], 
                       by=list( population=wings$population, 
                                sex=wings$sex, 
                                line=wings$line),
                       FUN=mean )

procoords <- wings_id[,4:99] * matrix( rep( rep( c(-1, 1), 48), 94), ncol= 96, byrow = T)

meanshape <- colMeans( procoords )

ZambiaEthiopia_diff <- colMeans( procoords[ wings_id$population == "Zambia", ] ) - 
  colMeans( procoords[ wings_id$population == "Ethiopia", ] )

pdf(file ="../outputs/Figure1C_MeanShapeDiffHighAltLowAlt.pdf", width = 12, height = 12) 
c_1 <- WingEffect( meanshape, ZambiaEthiopia_diff , 
                 ZambiaEthiopia_diff,
                 wingcol=c("black","black", "gray60" ),
                 scale.factor = 2,
                 scale.display = FALSE,
                 wingframe = FALSE,
                 winglwd=c(5, 5, 5))

dev.off()


# loading cropped version of wingplot as a pdf
c_2 <- image_read_pdf("../outputs/Figure1C_CroppedMeanShapeDiffHighAltLowAlt.pdf")

# turing wingplot pdf into ggplot friendly object
c <- ggdraw() + draw_image(c_2, scale = 1)

###############################################################


################# Wing shape raw data - PCA ###################

wing_PCs <- prcomp(wings[,10:105])
wingsPC <- data.frame(wings, wing_PCs$x[,1:58])

# PC1 vs PC2 - Plot for check only, not included in paper
d <- ggplot(data = wingsPC, aes(x = PC1, y = PC2, color = population, shape = population))+
  geom_point(show.legend = TRUE, size = 8)+
  theme(legend.position = "bottom", legend.box = "vertical") +
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude")) +
  scale_shape_discrete (name  ="Population", 
                        labels = c("High-Altitude", 
                                   "Low-Altitude"))

d

# PC1 vs Size - Plot for check only, not included in paper
e <- ggplot(data=wingsPC, aes(x=CS_rescaled, y = PC1, color= population, shape = population))+
  geom_point(show.legend = TRUE, size = 8) +
  xlab("Wing Size")+
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        legend.background = element_rect(fill = alpha("white", 0), size=0.5))+
  scale_colour_manual(values = c("black", "gray60"),
                      name  ="Population", 
                      labels = c("High-Altitude", 
                                 "Low-Altitude")) +
  scale_shape_discrete (name  ="Population", 
                        labels = c("High-Altitude", 
                                   "Low-Altitude"))
e

###############################################################


#### Microenvironmental canalization - Wing Shape analysis #####

# converting data into geomorph-friendly format
wings_land<- wings[,10:105]
wings_descr <- wings[,c(3, 4, 5, 107)]
wings_csize <- wings$CS_rescaled

#converting landmarks into 3D array because geomorph function uses this type of input
wings_land_3D <- arrayspecs(wings_land, 48, 2)

#converting to list
wings_list_3Darray <- list(wings_descr = wings_descr, 
                           wings_land_3D = wings_land_3D, 
                           wings_csize = wings_csize)

#converting to geomorph data frame
gdf <- geomorph.data.frame(shape = wings_land_3D,
                           population = wings_descr$population, 
                           sex = wings_descr$sex,
                           line = wings_descr$line,
                           size = wings_csize)


## Wing Shape Model: using procD.lm from geomorph library
wings_allomAnova <- procD.lm(shape ~ size * sex * (population/line), 
                             data = gdf, iter = 5000,
                             RRPP=TRUE, print.progress = FALSE)

summary(wings_allomAnova)

# updating the residuals to account for nested design
wings_allomAnova_nested <- anova(wings_allomAnova, 
                                 error = c("Residuals", "Residuals",
                                           "population:line", "Residuals",
                                           "Residuals", "Residuals",
                                           "Residuals","Residuals", 
                                           "Residuals","Residuals",
                                           "Residuals")) 

# Table S5 - Wing shape Manova
summary(wings_allomAnova_nested)

###############################################################


################ Figure 1 #################
bottom_row <- plot_grid(b, c, labels = c('B', 'C'), label_size = 65)

Figure1_test <- plot_grid(a, bottom_row , ncol = 1, labels = c('A',''), label_size = 65, rel_heights = c(2, 1))

pdf(file = "../outputs/PaperFigures/Figure1_WingSize_WingShape_Landmarks.pdf", height = 32, width = 36)
Figure1_test
dev.off()
##############################################