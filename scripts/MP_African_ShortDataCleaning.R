####################################################################################
################## Figure 1, 2, 3 ORIGINAL DATA CLEAN UP SCRIPT ####################
####################################################################################

## This script is to check and if necessary clean 
##(or remove) data based on notes in readme_specimeninfo.txt. 
##Mostly issues with possible mis-labeling for a few samples.


############ Loading Libraries ###########

library(data.table)

############ Read in the data ###########

wings <- read.csv( "../data/MP_African2016_AllPoints_Output_Landmarks+Veins+Outline.csv", 
                   header=TRUE )

############# Create population variable ############
##Creating a population level variable for Ethiopia and Zambia
## First: extract first character only
pop_dummy <- substr(wings$line, 1,1)
wings$population <- ifelse(pop_dummy == "z", "Zambia", "Ethiopia")
wings$population <- factor(wings$population)

############# Fix the size scaling because of recropping ############
# We (MP and ID) Calculated how the cropping and resizing effects the scale.
# For "Wings" scale is based on the resized to 632x480 images.

# For original (raw) images at 1360x1024, once resized to 632x480, 144px/mm or 0.007 mm/px

# Batch 1-3 was cropped to 1087x703. Once resized to 632x480, we end up with 182 px/mm or (for software) 0.0055 mm/px

# Batch 4 was cropped to 1046x788. Once resized to 632x480 we end up with 188 px/mm or 0.0053 mm/px

# Batch 5 and most of 6 was cropped to 907x683. Once resized to 632x480 we end up with 218 px/mm or 0.0046 mm/px

# A few batch 6 images had to be cropped to 950x726.Once resized to 632x480 we end up with 206 px/mm or 0.0049 mm/px. These are African_M_EF131N_014, African_M_EF131N_015, African_M_EF131N_016, African_M_EF131N_017 and African_M_EF131N_018. These are called batch 7.

# So let's make a vector to do centroid size corrections based on this.

# First create a vector for SplineBatch
SplineBatch <- unique(wings$SplineBatch)

scaleFactor <- c(182, 182, 182, 188, 218, 218, 206)

dummy_variables <- data.frame(SplineBatch, scaleFactor)

wings <- merge(wings, dummy_variables, 
    by.x= "SplineBatch", by.y = "SplineBatch")
    
#And now recompute our centroid size by dividing by scale factor and multipling by 144
wings$CS_rescaled <- with(wings, (CS/scaleFactor)*144)


############# Fix the lines with unclear names ############
##There was some issues with names of lines. Here we fix those issues.
# The ef39n_95n females are ef39n females. rename in full dataset. 

fix1 <- ifelse(wings$line =="ef39n_95n", "ef39n", as.character(wings$line))
wings$line <- factor(fix1)

## Images with the "-1.tif" added to "ef96n males" 
## are actually are 86n males,not ef96n males
check1<-wings$FileName[(grepl("ef96", wings$FileName)) & (grepl("-1.tif", wings$FileName))]
index1 <- (grepl("ef96", wings$FileName)) & (grepl("-1.tif", wings$FileName))
fix2 <- ifelse(index1, "86n", as.character(wings$line))
wings$line <- factor(fix2)

############ Remove lines with small sample sizes ############

## Remove lines with very small sample sizes.
## This includes lines with very few individuals
## and lines with only individuals from one sex
mydf <- wings

df2 <- mydf[with(mydf, as.logical(ave(as.character(line), line, 
                               FUN = function(x) length(x) > 25))),]
dim(mydf)
dim(df2)
df2$line <- droplevels(df2$line)
df2$sex <- droplevels(df2$sex)

############ Add "ef" to lines that did not have it ###############
## Some of the lines did not have an ef in front of the line number.
## These lines only have numbers as line names. 
## Here we add "ef" to these lines.
lines_without_ef <- grepl(pattern = "^[0-9]", df2$line)
lines_without_ef2 <- ifelse(lines_without_ef, 
    paste("ef",df2$line, sep = ""),
    as.character(df2$line))
lines_without_ef2 <- as.factor(lines_without_ef2)
df2$line <- lines_without_ef2

wings <- df2


################# Remove ef134n because they seem to be size outliers #################

wings <- subset(wings, subset = line != "ef134n", drop = TRUE)

wings$line <- droplevels(wings$line)


################### Data Clean up is Done ########################## 