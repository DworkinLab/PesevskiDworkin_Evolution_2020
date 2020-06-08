####################################################################################
#################### Figure 4 CELL COUNTS DATA CLEAN UP SCRIPT #####################
####################################################################################

######## Load Libraries ##########
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)

###################### Load Data #######################
cells_raw <- read.csv("../data/CellCounts_AfricanExperiment.csv", h=T)

##################### Clean Up Data #######################
### Creating and Fixing Descriptor variables

# change the names of the variables
cells_names <- c("wing_ID", "wing_area","area_location", "trichome_num", "trichome_den")
names(cells_raw) <- cells_names

# Change area_location to a factor
cells_raw$area_location <- factor(cells_raw$area_location) 

#Make a new factor called wing_regions
levels(cells_raw$wing_area) <- c("A","A","B","B","B", "C","C", "C", "E", "E", "F", "F","F","G", "G", "G")

# creating a variable labelling each individual measurement by a letter and number
cells_raw$wing_regions <- with(cells_raw, 
                           interaction(area_location, wing_area,  sep ="", drop=T))

# Duplicating the wing_ID column so that we can separate the file name into descriptor variables
cells_raw$wing_ID2 <- cells_raw$wing_ID

# removing file extension first:
cells2 <- separate(cells_raw, wing_ID2, c("wing_ID3", "extension"), sep = "\\.")


# creating predictor varaibles from file name:
cells3 <- separate(cells2, wing_ID3, c("initials", "sex", "line1", "line2", "individual"), sep="\\_") # separating file name

cells4 <- cells3[,-c(7,12)] # removing irrelevant columns

cells4$population <- cells4$line1 # creating a population variable

cells4[,c(7:11)] <- lapply(cells4[,c(7:11)], factor) # making new predictor variables into factors

levels(cells4$population) <- c("Ethiopia", "Ethiopia", "Zambia") # renaming population levels

levels(cells4$line1) <- c("EF", "EF", "Z") # fixing line1 levels so that they are all upper case

cells4$line <- interaction(cells4[,8], cells4[,9], drop=TRUE, sep ="") # combining line1 and line2 to create a single line name

cells5 <- cells4[,-c(8,9)] # removing irrelevant columns

cells <- cells5 # simplifying data frame name.


### Converting cell density within 75x75 px area to cells per mm^2

# converting the cell counts into cells per mm^2
cells$Area_mm2 <- rep(0.0065, 17872 )

cells$Dens_per_mm2 <- cells$trichome_num/cells$Area_mm2

cells$cell_area <- 1/cells$Dens_per_mm2

### Removing lines that we know are mislabelled
# removing E134 because they are size and cell density outliers
cells <- subset(cells, subset = line != "EF134N", drop = TRUE)

cells$line <- droplevels(cells$line)

