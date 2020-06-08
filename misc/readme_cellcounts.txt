This is a readme file for the Cell Count African II project

This study uses a re-imaged subset of the African II project wings at a much higher resolution, for trichome counts using Fijiwings. 

The subset includes 10 random lines from each Ethiopia and Zambia chosen using the sample() function in R. These lines are: "ef7n"   "ef15n" "ef122n" "ef1n"   "ef136n" "73n"    "31n"    "ef25n"  "ef54n" "ef117n" "z217n" "z216n" "z251n" "z332n" "z360n" "z357n" "z461n" "z507n" "z186n" "z403n”.
**** ef7n slide was damaged during the re-imaging process so a line ef134n was chosen to replace it

Another subset of the 5 most and least variable (in size) lines was chosen as well. This list includes: 

bottom 5: z383n (0.04874949); z322n (0.05034296); z507n (0.05052992); z254n (0.05054911); z311n (0.05129699)

top 5: ef131n (0.07399871 ); ef98n (0.07048139); 16n (0.06965806); ef119n (0.06918701); ef19n (0.06903155)


A total of 1091 images were imaged using the Olympus BX43 system microscope with Olympus DP80 Camera, using the cellSens Standard (version 1.14) imaging software. 
The images were taken using the 4x objective (FN26.5) for a Total magnification of 40x. The images were taken at 4080 × 3072 resolution at 8 bit depth in RBG colour.

**********************************
***Repeitability experiment:***
One Ethiopia (ef15n) line and one Zambia (z217n) line were used to test the repeatability of trichiome (cell) count software FijiWings (version 2.2). Three observations were taken per each region of the wing. The regions are:
A - Between L1 and L2
B - Between L2 and L3
C - Anterior side of Posterior cross vein between L3 and L4 vein
D - Posterior side of Posterior cross vein between L3 and L4 vein
E - Posterior side of anterior cross vein between L4 and L5
F - Anterior side of anterior cross vein between L4 and L5
G - Between L5 and the edge of the wing

The observations were taken at random in the centre of each region except for region D. A 75px square (total area 22.5 square kilopixels) was used for the observation. The data was summarized in an Excel spreadsheet. The data was analyzed using R (version 3.2.3) to test the variability between each observation and to compare it to the variability between line, population, region and sex. 

***The variability between observations was almost 0 so we concluded that repeatability is very high. 

*********************************

***Cell count experiment***

Trichome number observations for the images for the above-mentioned lines (10 random and 5 most and least variable) were taken using FijiWings (version 2.2). The observations were taken in different locations of the different wing areas. The same wing area designation is used as above. 
For area A two observations were taken: 1 is the proximal and 2 is the distal observation
For area B three observations were taken: 1 is the proximal and 2 is the central and 3 is the distal observation
For area C three observations were taken: 1 is the proximal and 2 is the central and 3 is the distal observation
For area E two observations were taken: 1 is the proximal and 2 is the distal observation
For area F three observations were taken: 1 is the proximal and 2 is the central and 3 is the distal observation
For area G three observations were taken: 1 is the proximal and 2 is the central and 3 is the distal observation

A 75px square (total area 22.5 square kilopixels) was used for the observation. The data was summarized in an Excel spreadsheet. The data was analyzed using R (version 3.2.3). 