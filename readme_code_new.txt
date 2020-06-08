This is a readme file for the environmental canalization project using populations of Drosophila melanogaster from low and high altitude populations from Africa

Description of how everything is organized:

Main Folders:

"data" - The folder named "data" contains all of the data files (.csv file extensions and .dat file extensions), with all of the relevant data for this project. These include:
	- African_Condition_2_BinaryFormatMutants_macrocanalization.csv - data file containing frequencies of mutation defects for flies raised at different temperatures
	- African_Mutants.xlsx - this is a data file containing the wings that contain mutation defects classified by the type of mutation that was observed
	- African_Mutation_frequencies.csv (and .xlsx) - data file containing frequencies of mutation defects for microcanalization experiment
	- CellCounts_AfricanExperiment.csv - data file containing cell density data
	- MP_African2016_AllPoints_Output_Landmarks+Veins+Outline.csv (and .dat) - data file containing wing size and shape data for the original microcanalization experiment
	- MP_AfricanCondition2017_EandZ_Output_Landmarks+Veins+Outline.csv (and .dat) - data file containing wing size and shape data for the macro canalization experiment (temperature manipulation)
	- MP_FA_only_AfricanCondition_ME_Output_Landmarks+Veins+Outline.csv (and .dat) - data file containing wing size and shape data for the fluctuating asymmetry experiment.



"misc" - The folder named "misc" contains the following files:
	- additional readme files (readme_specimeninfo.txt andreadme_cellcounts.txt) describing how the original data set(microenvironmental canalization) and the cell counts data set were collected along with notes on any issues with labelling, scaling etc. 
	- In this folder there is the file "Population_info" that describes the collection record for the strains from each population. 
	- There are two image files that represent the observations regions where the cell count measurements were taken (Observation_positionsandregions.png, Observation_positionsandregions_cropped.png). 
	- There is a pdf file of a wing image representing the landmarks and semi-landmark positions used in the shape analysis. 
	- There are two code files with custom functions created by our lab, one used to produce wing shape plots (WINGPLOTSOURCE.R) and another one used to calculate Levene's deviates (ID_LeveneStat_V1_2016).




"outputs" - The folder named "outputs" contains the individual plots that make up the panels of main and supplementary figures that are outputted from the analysis code files. Within this folder there is a folder called "PaperFigures" that contains the full panelled figures that are included in the manuscript. 




"scripts" - The folder named "scripts" contains the R code files that are used to clean up the data and perform analysis. Below is a full description of each script file.
-------------------------------------------------------------------------------------------------------------------------------


Description of R code files in the scripts folder:

Clean up code files:

1. MP_African_ShortDataCleaning.R - the code file cleaning the original data set (microenvironmental canalization).

- this code requires the following files:
	* Data csv file: MP_African2016_AllPoints_Output_Landmarks+Veins+Outline.csv - in "data" folder

- this code is called in the following scripts:
	* MP_AfricanCanalization_Figure1_WingSize_Shape.R
	* MP_AfricanCanalization_Figure2_WingSize_CV_LevDev_WingDefects.R
	* MP_AfricanCanalization_Figure3_WingShape_TotVar_Ecc_Defects.R

- this code is separated into 6 main sections
	1. Create population variable
	2. Fix the size scaling because of recropping - because the data was phenotyped in multiple batches, some of the batches had to be recropped in order maximize the number of wings that are splined. The scale that was originally entered was not changed to correspond to the recropped images. We therefore had to fix the scaling issues by multiplying with the appropriate scaling factors in the code. 
	3. Fix the lines with unclear names - when imaged, some of the wings were mislabelled and these issues are resolved in this section (notes for the appropriate names are in file readme_specimeninfo.txt)
	4. Remove lines with small sample sizes - some lines had very few wings so these were removed
	5. Add "ef" to lines that did not have it - when imaged the population designations for some of the high-altitude Ethiopian wings were missed. These are fixed in this section
	6. Remove ef134n because they seem to be size outliers - this line was removed because it was confirmed that the genotype was mislabelled.

2. MP_African2016_CellCounts_DataCleanup.R - the code file cleaning the cell counts data set.
- this code requires the following files:
	* Data csv file: CellCounts_AfricanExperiment.csv - in "data" folder

- This code is called in the following script:
	* MP_AfricanCanalization_Figure4_CellDensity_andCV.R

- This code is separated into 3 main sections:
	1. Creating and fixing descriptor variables from the file names (population, sex, line etc.)
	2. Converting the cell density from number of cells per 75 x 75 px to cells per mm^2
	3. Removing a strain that was phenotyped that was later confirmed it was mislabelled.

3. MP_MacroEnvironmental_Canalization_Cleanup.R - the code file cleaning the second data set with temperature manipulation (macroenvironmental canalization)

- this code requires the following files:
	* Data csv file: MP_AfricanCondition2017_EandZ_Output_Landmarks+Veins+Outline.csv - in "data" folder

- This code is called in the following script:
	* MP_AfricanCanalization_Figure5_Macroenvironmental_Canalization.R

- This code is separated into 2 main sections
	1. Creating and fixing descriptor variables from file names (population, sex, line, temperature etc.)
	2. Removing data that is either from mislabelled lines, lines with low sample sizes, or data that is irrelevant because it was from flies raised at 15% and 5% nutritional dilutions (for another experiment)

4. MP_African_Canalization_FA_cleanup_outlier.R - the code file cleaning the third data set with left and right wings (fluctuating asymmetry)

- this code requires the following files:
	* MP_FA_only_AfricanCondition_ME_Output_Landmarks+Veins+Outline.csv

- This code is called in the following scripts:
	* Data csv file: MP_AfricanCanalization_Figure6_FA.R - in "data" folder

- This code is separated into two main sections:
	1. Fixing the descriptor variables to have consistent labels for temperature, adding numeric temperature variable, adding a population variable, and subsetting the wings that were measured twice for the measurement error estimation.
	2. Calculating the different FA measures for wing size (FA1 and FA8) and shape (PD_Left_Right). This section includes a hashed out outlier check section for the FA1 and FA8 measures that was used in order to remove outliers values that were 
greater than 3x the SD of each measure. 

-------------------------------------------------------------------------------------------------------------------------------

Analysis script files:

1. MP_AfricanCanalization_Figure1_WingSize_Shape.R - This code file includes the analysis looking at the differences in wing size and shape between the different populations and sexes and code to produce Figure 1.

- This code requires the following files
	* Data cleaning script: MP_African_ShortDataCleaning.R - in the "scripts" folder
	* Local functions: WINGPLOTSOURCE.R - in the "misc" folder
	* Image file: Figure1B_Landmarks.png - in the "outputs" folder
	* Image file: Figure1C_CroppedMeanShapeDiffHighAltLowAlt.pdf - in the "outputs" folder

- This code is separated into 7 main sections
	1. Code creating Figure 1A using the raw data
	2. Linear mixed model testing the effect of population and sex on wing size
	3. Loading the image for figure 1B
	4. Creating the wing plot for figure 1C
	5. PCA analysis of wing shape - not included into the manuscript, for visualizing shape variation only
	6. Multivariate Linear model testing the effect of population and sex on wing shape
	7. Code producing the panelled Figure 1.


2. MP_AfricanCanalization_Figure2_WingSize_CV_LevDev_WingDefects.R - This code file includes the analysis comparing the different within-line measures of variability (microenvironmental variation) between the two population for wing size and looking at associations between the different measures of variability with wing defects. This code file is used to produce Figure 2, Supplementary Figure 1, Supplementary Figure 3 and Supplementary Figure 4.  

- This code requires the following files:
	* Data cleaning script: MP_African_ShortDataCleaning.R - in the "scripts" folder
	* Local functions: ID_LeveneStat_V1_2016.R - in the "misc" folder
	* Mutant frequency data file: African_Mutation_frequencies.csv - in the "data" folder
	* Image file showing examples of wing defects: Example_WingDefects.png - in the "outputs" folder

- This code is separated into 11 main sections
	1. Code calculating the within-line CV for each line, a generalized linear model testing the effect of population on within-line CV
	2. Code calculating within-line Levene's deviates, a generalized linear model testing the effect of population on within-line Levene's deviates
	3. Generalized model testing the effect of population and sex on mutation frequencies.
	4. Code examining at the associations between within-line CV and within-line mutation frequencies 
	5. Code examining at the associations between within-line Levene's deviates and within-line mutation frequencies 
	6. Code examining the relationship between CV and Levene's deviates
	7. Code calculating the Coefficient of Genetic Variation to look at the difference between the two populations. This is not included in the paper, it is only for checking.
	8. Code producing Figure 2
	9. Code producing Supplementary Figure 1
	10. Code producing Supplementary Figure 3
	11. Code producing Supplementary Figure 4


3. MP_AfricanCanalization_Figure3_WingShape_TotVar_Ecc_Defects.R - This code file includes the analysis comparing the different within-line measures of variability (microenvironmental variation) between the two population for wing shape and looking at associations between the different measures of variability for wing shape with wing defect frequency. This code file is used to produce Figure 3, Supplementary Figure 2, Supplementary Figure 5 and Supplementary Figure 6.

- This code requires the following files:
	* Data cleaning script: MP_African_ShortDataCleaning.R - in the "scripts" folder
	* Mutant frequency data file: African_Mutation_frequencies.csv - in the "data" folder

- This code is separated into 8 main sections
	1. Code calculating the within-line total variance, eccentricity, rSDE and rSDE2 for each line.
	2. Generalized linear model testing the effect of population on within-line total variance, eccentricity, rSDE and rSDE2
	3. Code examining associations between within-line total variance, eccentricity, rSDE and rSDE2 with within-line mutation frequencies.
	4. Code examining the morphological disparity within and between the two populations. 
	5. Code producing Figure 3
	6. Code producing Supplementary Figure 2
	7. Code producing Supplementary Figure 5
	8. Code producing Supplementary Figure 6


4. MP_AfricanCanalization_Figure4_CellDensity_andCV.R

- This code requires the following files
	* Wing size and shape data cleaning file: MP_African_ShortDataCleaning.R - in "scripts" folder
	* Cell density data cleaning file: MP_African2016_CellCounts_DataCleanup.R - in "scripts" folder
	* Image file outlining cell count measurement locations: image18.png - in "misc" folder
	* Analysis code file: MP_AfricanCanalization_Figure2_WingSize_CV_LevDev_WingDefects.R - in "scripts" folder
	* Analysis code file: MP_AfricanCanalization_Figure3_WingShape_TotVar_Ecc_Defects.R - in "scripts" folder

- This code is separated in 10 main sections
	1. Making Figure 4A image into an R friendly object
	2. Linear mixed models testing the effect of population and sex, (wing regions/areas in some models) on wing cell density
	3. Code calculating the within-line CV for wing cell density calculated two ways: first by calculating the individual CV and taking the mean of each line and second by calculating the line CV. 
	4. Code calculating the within-line Levene's deviates for wing cell density
	5. Code examining the relationship between CV and Levene's deviates for wing cell density
	6. Code examining the relationship between cell density CV and wing size CV, wing shape total variance, eccentricity
	7. Code examining the association between cell density CV and wing defects frequency
	8. Code producing Figure 4
	9. Code producing Supplementary Figure 7
	10. Code producing Supplementary Figure 8


5. MP_AfricanCanalization_Figure5_Macroenvironmental_Canalization.R

- This code requires the following files
	* Data clean file: MP_MacroEnvironmental_Canalization_Cleanup.R - in "scripts" folder
	* Data csv file: African_Condition_BinaryFormatMutants.csv - in "data" folder
	* Local functions: WINGPLOTSOURCE.R - in the "misc" folder

- This code is separated into 16 main sections:
	1. Getting the average wing size and shape for each population for plotting purposes
	2. Removing any lines raised at specific temperatures that have N less than 20 (both males and females)
	3. Linear mixed model testing the effect of population, sex and temperature on wing size. Looking at the wing size reaction norms for each population and sex due to changes in developmental temperature.
	4. Linear model testing the effect of size, population, sex and temperature on wing shape. 
	5. PCA of the wing shape data
	6. Code cleaning the mutant frequencies in wings raised at multiple temperatures and running a generalized linear mixed model testing the effect of temperature, population and sex on frequency of wing defects.
	7. Calculating the within-line CV and Levene's Deviates for wings raised at different temperatures and linear mixed models testing the effect of population and temperature on within-line CV and Levene's Deviates
	8. Interlude code for checking correlations of within-line CV and Levene's deviates between lines raised at different temperatures - for internal checks only this is not included in the manuscript
	9. Calculating the within-line total variance, eccentricity, rSDE and rSDE2 for wing shape for wings at different temperatures. Running generalized linear mixed models testing the effects of temperature and population on within-line measures of variability for wing size
	10. Interlude code for checking correlations of within-line total variance and eccentricity for wing shape between lines raised at different temperatures - for internal checks only this is not included in the manuscript
	11. Code examining the association between within-line CV and Levene's deviates for wing size and within-line total variance, eccentricity, rSDE and rSDE2 for wing shape with within-line frequency of wing defects.
	12. Creating wing shape plots to compare mean shape differences between high and low altitude wings at different temperatures. Calculating the PD (difference in magnitude) between the high and low altitude mean shapes at different temperatures and calculating the correlation (difference in direction) between the high and low altitude mean shape vectors at different temperatures
	13. Code producing Figure 5
	14. Code producing Supplementary Figure 9
	15. Code producing Supplementary Figure 10
	16. Code producing Supplementary Figure 11


6. MP_AfricanCanalization_Figure6_FA.R

- This code requires the following files:
	* Data cleanup file: MP_African_Canalization_FA_cleanup_outlier.R - in "scripts" folder
	* Analysis code file: MP_AfricanCanalization_Figure5_Macroenvironmental_Canalization.R - in "scripts" folder

- This code is separated into 9 main sections:
	1. Generalized linear mixed models testing the effects of temperature sex and population for FA1, FA8 and PD_Left_Right. 
	2. Estimating Measurement error for FA for wing size and wing shape
	3. Testing differences in wing shape FA using Geomorph function bilat.symmetry where we first remove the directional asymmetry and look at the differences in fluctuating asymmetry between the populations at different temperatures
	4. Redoing the analysis using bilat.symmetry separately for each population in order to avoid assuming similar DA - for checking only, not included in manuscript.
	5. Code examining the associations between measures of FA for wing size and shape with within-line measures of variation for wing size and shape at different temperatures.
	6. Code producing Figure 6
	7. Code producing Supplementary Figure 12
	8. Code producing Supplementary Figure 13
	9. Code producing Supplementary Figure 14
