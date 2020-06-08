####################################################################################
#################### Figure 5 TEMPERATURE DATA CLEAN UP SCRIPT #####################
####################################################################################


############### Loading Libraries ##############
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(MuMIn)
#################################################

pref_theme <- theme_classic() + 
  theme(text = element_text(size = 28), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
theme_set(pref_theme)

################## Reading Data #################
wings_raw <- read.csv("~/Dropbox/MariaPesevski/AfricanCondition2016/data/MP_AfricanCondition2017_EandZ_Output_Landmarks+Veins+Outline.csv")

#################################################

############### Cleaning Data ###################

### Creating descriptor varaibles from file name.

#Creating a separate file column for extraction of information from file name 
wings_raw$File2 <- wings_raw$File

#Extracting the factors from the file names
wings2 <- separate(data = wings_raw, col = File2, into= c("file1", "file2", "file3", "file4", "file5","wingID"), sep = "\\\\")

wings3 <- wings2[,-c(1,3:6,12,110:114)] # removing irrelevant columns

wings4 <- separate(data = wings3, col = wingID, into= c("initials", "line", "condition", "replicate", "sex","individual"), sep = "_")

wings_clean <- wings4

wings_clean[,104:109] <- lapply(wings_clean[,104:109], factor) # making descriptor variable columns factors

#creating the nutrition and temperature factors from the condition factor. 
nutrition <- wings_clean$condition

levels(nutrition) <- c("f100", "f15", "f5", "f100", "f100", "f100")

temperature <- wings_clean$condition

levels(temperature) <- c("t24", "t24", "t24","t18", "t24", "t28")

wings_clean <- cbind(wings_clean, nutrition, temperature)

#Creating a numeric temperature for analysis
new_temp <- wings_clean$temperature
levels(new_temp) <- c("24", "18", "28")
temp_num <- as.numeric(levels(new_temp))[new_temp]
wings_clean <- cbind(wings_clean, temp_num)
wings_clean$temp_zeroed <- wings_clean$temp_num - 18


wings <- wings_clean

# cleaning up the line levels
levels(wings$line) <- c("e112",  "e117",  "e119",  "e122",  "e131",  
                        "e134", "e136", "e15",   "e19",   "e39",   
                        "e54",   "e65", "e73",   "e119", "e16",  
                        "e217", "e98",  "z124","z159",  "z160",  
                        "z186",  "z216",  "z217",  "z251", "z254",  
                        "z311",  "z360",  "z366",  "z367",  "z455", 
                        "z461")

# creating a population factor

pop_dummy <- substr(wings$line, 1, 1)

wings$population <- ifelse(pop_dummy == "z", "Zambia", "Ethiopia")

wings$population <- factor(wings$population)

# releveling the temperature factor
wings$temperature <- relevel(wings$temperature, "t18")

with(wings, table(line:sex, condition))


### Getting rid of data that is irrelevant

#getting rid of observation that only had observation for a single sex
missing_index <- with(wings, line:nutrition %in% c("e122:f15", "e136:f15"))
wings <- wings[(!missing_index), ]
with(wings, table(line:sex, condition))
dim(wings)

#getting rid of e134 because I believe it is mislabelled (by the Pool lab)
wings <- subset(wings, subset = line != "e134", drop = T)

#getting rid of e217 and z217 because its likely they got mixed up during egg laying.
wings <- subset(wings, subset = line != "e217", drop = T)
wings <- subset(wings, subset = line != "z217", drop = T)
wings$line <- droplevels(wings$line)

# removing lines with no observation for certain sexes and temperatures
wings <- subset(wings, subset = line != "e54", drop = T)

wings <- subset(wings, subset = line != "z160", drop = T)

wings <- subset(wings, subset = line != "z461", drop = T)

wings$line <- droplevels(wings$line)

# Getting rid of data for flies that were raised at 15% and 5% food because it is not relevant.
wings_temp <- as.data.frame(subset(wings, nutrition == "f100", drop = T))
