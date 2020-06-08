This is a readme file for the African II project.

This is a follow-up study to Pitchers et al 2013.

Unlike the original study, many of these lines are more recently collected. These were all reared in the Pool lab (by Justin Lack) at the University of Wisconsin. They were sent in Ethanol (70%) to the Dworkin lab at MSU.

The wings were dissected by Megan Cermak (slides 1-?) and Christian Marier.

The wings were imaged by Ian Dworkin (April 2015). They were imaged at MSU on an Olympus BX51 microscope (with 2X objective for a total of 20X magnification). We used the Olympus DP controller software. I used the "Dworkin_NewWingProfile2014_2X.env" environment file for conditions (basically vanilla settings for digital contrast and sharpness settings". This also leaves a 2mm scale bar on the bottom right of each image.

This also put the contrast/balance meter (SPOT) and the fine focus (FF) grids on. I put the corner of the SPOT meter (bottom left 15%) over the L1-L2 intersection. For FF, I check it for each slide at the L5-PCV intersection, as well as L5-margin.

Image Naming convention

Project_Sex_Line_Individual

i.e.
African_F_EF135_01

Note for imaging. I put the auto-balance at the intersection of L1-L2.

Note: For EF131N (Dissected by Megan) both slides are labeled as male.However based on position in slide  box (Megan always had it F then M) and size of the wings. This issue was sorted out. 

Note: There are two slides labeled EF95N Males. One was placed right after EF39N females. As there is no corresponding slide for EF39N males, it is likely that this is the actually EF39N males. I have labeled images for it as African_M_EF39N_95N as a precaution.

I have maintained the EF95N males slide that I think is correct in the normal fashion (African_M_Ef95N_).

ONce splined, fit a CVA to both lines (and sexes, to double check).


Note: EF32N Males. Does not seem to have a counterpart female slide.


Note: On April 22nd 2015 (Starting with EF_M_73N) I changed one aspect of the "environmental settings file. Instead of a large 2mm scale bar, I switched to a smaller 500uM scale bar, and put it further down in the bottom right corner. I renamed the environment file as "Dworkin_NewWingProfile2014_2x_SmallerScale.env".


Note: Right after (EF)86N_F there was a slide of EF 96N males. Possible that these were EF 86N males (?). So far there was not another EF 96N males slide (but there was a female one). Use a discrimant function with 86N and 96N to sort this out.

Note: I also noticed that I (ID) forgot to put "EF" in front of a few of the lineage identifiers. However, so far all lines have begun with this, so it will not lead to any problems with consistency. Latser lines do have different identifiers, so care will be taken!

Note: African EF 132N. Only 1 male and 2 females. Sizes of both sexes are really similar.Double check using CVA that there is no evidence of reversal of sexual dimorphism).

NOTE: Ask Christian Marier about the slide notation. i.e. 2:397N  (Which I am writing down as 2-397N). Based on the notes from the pool lab it is not a "2" but a z (for Zambia). So I will change all of these later.

The wings were splined by Maria Pesevski (May 2016). They were splined at McMaster University and at Maria’s home. Wing Machine (Version 3.7) was used for splining and CPReader (Version ___) was used for superimposition. 


The images were resized to 632x480 for landmarking and 316x240 for splining using GIMP 2(Version 2.8). Landmarking was done using tpsDig (Version 2.16).



Batches Batch4_recropped, Batch5_Recropped and Batch6 were recropped because the wings would not spline. 

For Batch4_recropped dimensions 1046x788 were used for cropping, for Batch5_recropped dimensions 907x683 were used and for Batch 6 all images were cropped to 907x683 except for African_M_EF131N_014, African_M_EF131N_015, African_M_EF131N_016, African_M_EF131N_017 and African_M_EF131N_018 which were cropped to 950x726 (because using 907x683 cut off wing). 

These last 6 are called “batch 7”. The cropping was done in ImageJ using the Canvas Size tool. 

They were then resized to 632x480 for landmarking and 316x240 for splining using GIMP 2(Version 2.8).


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RESCALING PROCESS

Since the different batches were cropped at different sizes the following are the scaling factors that we used to rescale centroid size:

For original (raw) images at 1360x1024, once resized to 632x480, 144px/mm or 0.007 mm/px

Batch 1: cropped to 1087x703; scaling factor: 182 px/mm or 0.0055 mm/px
Batch 2: cropped to 1087x703; scaling factor: 182 px/mm or 0.0055 mm/px
Batch 3: cropped to 1087x703; scaling factor: 182 px/mm or 0.0055 mm/px
Batch 4: cropped to 1046x788. Scaling factor: 188 px/mm or 0.0053 mm/px
Batch 5: cropped to 907x683. Scaling factor: 218 px/mm or 0.0046 mm/px
Batch 6: cropped to 907x683. Scaling factor: 218 px/mm or 0.0046 mm/px
Batch 7: cropped to 950x726. Scaling factorL 206 px/mm or 0.0049 mm/px

The rescaling equation was as follows:

(CS/scaleFactor)*144

CS - centroid size
scaleFactor - the scaling factor appropriate for each batch (based on the cropping)
144 - is the scaling factor for uncropped images

We created a new variable called CS_rescaled with these new rescaled sizes.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



There was a total of 2182 raw images, out of those a total of 2052 were splined. The other 130 wings were rejected because of damage and the inability for Wing machine to spline them. 

Centroid size was calculated using all points, aligned points and unaligned points. 

NOTE: For the EF96N and the z322n lines there is two types of wing images with and without -1 after the numbering. These are different images with the same name except for that “-1” notation. 

Manual edits of the data frame prior to analysis:
NOTE: Some wings were splined twice (first time it incorrectly indicated these failed to spine without cropping). This generated duplicate entries. These duplicateswere checked and removed from the “MP_African2016_AllPoints_Output_Landmarks+Veins+Outline.csv” file.

NOTE: The 6 individual specimens noted above that were cropped differently (Batch 7) were manually edited to be called Batch 7 to indicate this in the “MP_African2016_AllPoints_Output_Landmarks+Veins+Outline.csv” file.


