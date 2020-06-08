####################################################################################
######################## Figure 6 FA DATA CLEAN UP SCRIPT ##########################
####################################################################################

########## Load Libraries and local functions #################
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
###############################################################


######################## Load Data ############################
#Size and shape data
wings_FA_raw <- read.csv("../data/MP_FA_only_AfricanCondition_ME_Output_Landmarks+Veins+Outline.csv")

###############################################################

#################### Set ggplot theme #########################

pref_theme <- theme_classic() + 
  theme(text = element_text(size = 65), 
        axis.title.y = element_text(margin = margin(t = 10, r = 20, b = 10, l = 10)),
        axis.title.x = element_text(margin = margin(t = 20, r = 10, b = 10, l = 10)),
        axis.text = element_text(margin = margin(t = 15, r = 15, b = 15, l = 15)))
theme_set(pref_theme)

###############################################################


####################### Data cleanup ##########################
### Fixing descriptor variablies

# Creating nutrition and temperature variables
nutrition <- wings_FA_raw$condition

levels(nutrition) <- c("f15", "f100", "f100", "f100")

temperature <- wings_FA_raw$condition

levels(temperature) <- c("t24", "t18", "t24","t28")

wings_FA_clean <- cbind(wings_FA_raw, nutrition, temperature)


#Creating a numeric temperature for analysis
new_temp <- wings_FA_clean$temperature
levels(new_temp) <- c("24", "18", "28")
temp_num <- as.numeric(levels(new_temp))[new_temp]
wings_FA_clean <- cbind(wings_FA_clean, temp_num)
wings_FA_clean$temp_zeroed <- wings_FA_clean$temp_num - 18

wings_FA <- wings_FA_clean

# creating a population factor

pop_dummy <- substr(wings_FA$line, 1, 1)

wings_FA$population <- ifelse(pop_dummy == "z", "Zambia", "Ethiopia")

wings_FA$population <- factor(wings_FA$population)

# releveling the temperature factor
wings_FA$temperature <- relevel(wings_FA$temperature, "t18")
wings_FA$population <- relevel(wings_FA$population, "Zambia")

# creating a separate data set with the wings phenotyped twice for measurement error calcualtions
wings_ME <- subset(wings_FA, subset = double == "d", drop = T)
wings_ME$inds <- interaction(wings_ME$line, 
                             wings_ME$condition,
                             wings_ME$replicate,
                             wings_ME$sex, 
                             wings_ME$individual)
###############################################################

####### Calculating FA measrues for wing size and shape #######

### 1. Making a data frame with size only were each individual has a wing size for L and R in separate columns
wings_SizeFA <- wings_FA[,c(2:8, 110:114, 109)]


wings_SizeFA_1 <- subset(wings_SizeFA, subset = double != "d", drop = T) # removing the wings that were phenotyped twice from original data set

wings_SizeFA_2 <- spread(wings_SizeFA_1, side, CS) # making L and R wing columns and putting the appropriate CS for each side. 

wings_SizeFA_3 <- na.omit(wings_SizeFA_2) # removing NAs

# creating a separate data set where CS is averaged between the two sides
wings_SizeFA_3.5 <- wings_SizeFA_3 

wings_SizeFA_3.5$meanCS <- rowMeans(wings_SizeFA_3.5[,12:13]) # taking the mean wing size for L and R wing for each specimen


wings_SizeFA_4 <- gather(wings_SizeFA_3, l, r, key = side, value = CS) # converting L and R wings abck into a single column turining it into long form

wings_SizeFA_4$side <- factor(wings_SizeFA_4$side)


### 2.Making a data frame where each individual has a wing shape and size for L and R in separate columns

wings_ShapeFA <- wings_FA[,c(2:8, 110:114, 13:109)]

wings_ShapeFA_1 <- subset(wings_ShapeFA, subset = double != "d", drop = T) # removing wings that were phenotype twice for ME

wings_ShapeFA_1$inds <- interaction(wings_ShapeFA_1$line, 
                                    wings_ShapeFA_1$condition,
                                    wings_ShapeFA_1$replicate,
                                    wings_ShapeFA_1$sex, 
                                    wings_ShapeFA_1$individual) # creating an individual column to label each specimen

wings_ShapeFA_2 <- wings_ShapeFA_1[order(wings_ShapeFA_1$inds),] # reordering the data frame the ind column

# merging the wing size and wing shape data frames 
wings_size_shape_FA <- merge(wings_SizeFA_4, wings_ShapeFA_1, by = c("line", "condition", "replicate","sex", "individual", "double","nutrition", "temperature", "temp_num", "temp_zeroed", "population", "side", "CS"))

### 3. Calculating the different FA measures for wing size

# calculating the difference in size between the right and left wing
wings_SizeFA_3.5$diff <- wings_SizeFA_3.5$r - wings_SizeFA_3.5$l
# calcilating the absolute value of the difference in size between left and right wing
wings_SizeFA_3.5$abs_diff <- abs(wings_SizeFA_3.5$r - wings_SizeFA_3.5$l)

# FA1 - calculating the mean difference in wing size between left and right wing for each line and sex
FA1_lines <- aggregate(wings_SizeFA_3.5["diff"],
                       c(wings_SizeFA_3.5["sex"],
                         wings_SizeFA_3.5["line"],
                         wings_SizeFA_3.5["nutrition"],
                         wings_SizeFA_3.5["temperature"],
                         wings_SizeFA_3.5["population"]),
                       mean)

# calculating the ln ratio for wing size between left and right wing
wings_SizeFA_3.5$ln_ratio <- log(wings_SizeFA_3.5$r/wings_SizeFA_3.5$l) 

# calculating the absoulte value of the ln ratio for wing size between the left and right wing
wings_SizeFA_3.5$abs_ln_ratio <- abs(log(wings_SizeFA_3.5$r/wings_SizeFA_3.5$l))

# FA8 - calculating the mean ln ratio in wing size between left and right wing for each line and sex
FA8_lines <- aggregate(wings_SizeFA_3.5[,16],
                       c(wings_SizeFA_3.5["sex"],
                         wings_SizeFA_3.5["line"],
                         wings_SizeFA_3.5["nutrition"],
                         wings_SizeFA_3.5["temperature"],
                         wings_SizeFA_3.5["population"]),
                       mean)
# Putting the two measures of FA into a single data frame
FA1_FA8 <- merge(FA1_lines, FA8_lines, by = c("sex", "line", "nutrition", "temperature", "population"))

#### Interlude: Checking for FA1 and FA8 outliers ####

# removing wings raised at 15% nurrition
wings_SizeFA_3.5_Temp <- subset(wings_SizeFA_3.5, nutrition %in% "f100")

#######################################################################
# # Checking outliers for FA1
# wings_SizeFA_3.5_Temp$FA1_diffs_OU <- abs(wings_SizeFA_3.5_Temp$abs_diff - mean(wings_SizeFA_3.5_Temp$abs_diff))
# 
# FA1_diffs_OU_Meds <- median(wings_SizeFA_3.5_Temp$FA1_diffs_OU)
# FA1_diffs_OU_Means <- mean(wings_SizeFA_3.5_Temp$FA1_diffs_OU)
# FA1_diffs_OU_SD <- sd(wings_SizeFA_3.5_Temp$FA1_diffs_OU)
# 
# FA1_diffs_OU_2timesSD <- 2*FA1_diffs_OU_SD
# FA1_diffs_OU_3timesSD <- 3*FA1_diffs_OU_SD
# FA1_diffs_OU_4timesSD <- 4*FA1_diffs_OU_SD
# 
# FA1meanPlus2SD <- FA1_diffs_OU_Means + FA1_diffs_OU_2timesSD
# FA1meanMinus2SD <- FA1_diffs_OU_Means - FA1_diffs_OU_2timesSD
# FA1meanPlus3SD <- FA1_diffs_OU_Means + FA1_diffs_OU_3timesSD
# FA1meanMinus3SD <- FA1_diffs_OU_Means - FA1_diffs_OU_3timesSD
# FA1meanPlus4SD <- FA1_diffs_OU_Means + FA1_diffs_OU_4timesSD
# FA1meanMinus4SD <- FA1_diffs_OU_Means - FA1_diffs_OU_4timesSD
# 
# 
# 
# bi <- ggplot(data = wings_SizeFA_3.5_Temp, aes(y=FA1_diffs_OU, x = line, colour = line, symbol = sex)) +
#   geom_hline(yintercept= c(FA1_diffs_OU_Meds,
#                            (FA1_diffs_OU_Means+FA1_diffs_OU_2timesSD), 
#                            (FA1_diffs_OU_Means-FA1_diffs_OU_2timesSD)), 
#              colour = "red") +
#   geom_point(size = 3)
# bi
# 
# bj <- qplot(data =wings_SizeFA_3.5_Temp, y = FA1_diffs_OU) +
#   geom_hline(yintercept= c(FA1_diffs_OU_Meds,
#                            (FA1_diffs_OU_Means+FA1_diffs_OU_2timesSD), 
#                            (FA1_diffs_OU_Means-FA1_diffs_OU_2timesSD)), 
#              colour = "red") +
#   geom_point(aes(y = FA1_diffs_OU, colour = FA1_diffs_OU > FA1meanPlus2SD), size = 3)+
#   scale_colour_discrete(name = "Outliers")
# 
# bj
# 
# # which(wings_SizeFA_3.5_Temp$FA1_diffs_OU > FA1meanPlus2SD)
# # which(wings_SizeFA_3.5_Temp$FA1_diffs_OU > FA1meanPlus3SD)
# 
# 
# # Checking Outlires for FA8
# wings_SizeFA_3.5_Temp$FA8_diffs_OU <- abs(wings_SizeFA_3.5_Temp$abs_ln_ratio - mean(wings_SizeFA_3.5_Temp$abs_ln_ratio))
# 
# FA8_diffs_OU_Meds <- median(wings_SizeFA_3.5_Temp$FA8_diffs_OU)
# FA8_diffs_OU_Means <- mean(wings_SizeFA_3.5_Temp$FA8_diffs_OU)
# FA8_diffs_OU_SD <- sd(wings_SizeFA_3.5_Temp$FA8_diffs_OU)
# 
# FA8_diffs_OU_2timesSD <- 2*FA8_diffs_OU_SD
# FA8_diffs_OU_3timesSD <- 3*FA8_diffs_OU_SD
# FA8_diffs_OU_4timesSD <- 4*FA8_diffs_OU_SD
# 
# FA8meanPlus2SD <- FA8_diffs_OU_Means + FA8_diffs_OU_2timesSD
# FA8meanMinus2SD <- FA8_diffs_OU_Means - FA8_diffs_OU_2timesSD
# FA8meanPlus3SD <- FA8_diffs_OU_Means + FA8_diffs_OU_3timesSD
# FA8meanMinus3SD <- FA8_diffs_OU_Means - FA8_diffs_OU_3timesSD
# FA8meanPlus4SD <- FA8_diffs_OU_Means + FA8_diffs_OU_4timesSD
# FA8meanMinus4SD <- FA8_diffs_OU_Means - FA8_diffs_OU_4timesSD
# 
# 
# bk <- ggplot(data = wings_SizeFA_3.5_Temp, aes(y=FA8_diffs_OU, x = line, colour = line, symbol = sex)) +
#   geom_hline(yintercept= c(FA8_diffs_OU_Meds,
#                            (FA8_diffs_OU_Means+FA8_diffs_OU_2timesSD), 
#                            (FA8_diffs_OU_Means-FA8_diffs_OU_2timesSD)), 
#              colour = "red") +
#   geom_point(size = 3)
# bk
# 
# bl <- qplot(data =wings_SizeFA_3.5_Temp, y = FA8_diffs_OU) +
#   geom_hline(yintercept= c(FA8_diffs_OU_Meds,
#                            (FA8_diffs_OU_Means+FA8_diffs_OU_2timesSD), 
#                            (FA8_diffs_OU_Means-FA8_diffs_OU_2timesSD)), 
#              colour = "red") +
#   geom_point(aes(y = FA8_diffs_OU, colour = FA8_diffs_OU > FA8meanPlus2SD), size = 3) +
#   scale_colour_discrete(name = "Outliers")
# bl
# 
# # which(wings_SizeFA_3.5_Temp$FA8_diffs_OU > FA8meanPlus2SD)
# # which(wings_SizeFA_3.5_Temp$FA8_diffs_OU > FA8meanPlus3SD)
#######################################################################

### Creating a new data frame with outliers removed based on analysis above that is hashed out
wings_SizeFA_NoOut_Temp <- wings_SizeFA_3.5_Temp[-c(13,  15,  16,  17,  37,  59, 106, 158, 159, 324, 417),]


### 4. Calculating FA measures for wing shape

# Procrustes distance between L and R for each specimen. 

left_data <- subset(wings_size_shape_FA, subset = side == "l", drop = T ) # extracting only left wings
right_data <- subset(wings_size_shape_FA, subset = side == "r", drop = T ) # extracting only right wings

# converting landmark data into vectors
left_vec <- data.matrix(left_data[,14:109])

right_vec <- data.matrix(right_data[,14:109])

# calculating the difference vector for the left and right wings
LR_Diff_Vec <- left_vec - right_vec

# taking the square root of the difference vector
LR_Diff_Vec_Sq <- (LR_Diff_Vec)^2

# adding the difference vector and the square of the difference vector along side the descriptor variables
LR_Diff_Vec_data <- cbind(left_data[,c(1:13, 110)], 
                          LR_Diff_Vec, LR_Diff_Vec_Sq )

# Calculating the Procrustes distance between left and right wings
vector_sums <- rep(NA, nrow(LR_Diff_Vec_Sq))

for(i in 1:nrow(LR_Diff_Vec_Sq)){
  vector_sums[i] <- sum(LR_Diff_Vec_Sq[i,])
}

PD_left_right <- rep(NA, length(vector_sums))
for(i in 1:length(vector_sums)){
  PD_left_right[i] <- sqrt(vector_sums[i])
}


# adding the PD_left_right column with the descriptor variables
PD_left_right <- cbind(left_data[,c(1:11, 13, 110)], PD_left_right)

###############################################################