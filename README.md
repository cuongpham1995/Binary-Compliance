# Overview
This is the codes for the simulations and data application of the binary compliance project. 

There are two parts in the simulations: 1) We know the true value of the sensitity parameter alpha; 2) we don't know the true value of the sensitivity parameter alpha. 

In each part, we consider four cases: 

Case 1: All nuisance parameters are correctly specified
Case 2: f(Z|X) are incorrectly specified
Case 3: f(A|Z,X) are incorrectly specified
Case 4: Q^z(X) are incorrectly specified

#Structure
There are four folders: Known Alpha, Unknown Alpha, Data Application and Heat Map Generator. 


##Known Alpha

The codes in the folder "Known Alpha" create the result in the section 4.1 in the paper. In this section, we know the true value of (alpha_0^-, alpha_0^+, alpha_1^-, alpha_1^+). 

"Data generator.R" simulates data set that is needed to run "All correct.R", "Faz wrong.R", "Fz wrong.R", "Q wrong.R". 
"All correct.R" corresponses to case 1. 
"Fz wrong.R" corresponds to case 2. 
"Faz wrong.R" corresponds to case 3. 
"Q wrong.R" corresponds to case 4. 
"Simulation known alpha.R" calls the functions from all the files above to run 500 simulations with different alpha values. The output of this file contains estimated correct classification rates, estimated value function (robust and non-robust), empirical value functions for four methods: proposed methods (IPW and MR); Eric Tchetgen's method, and OWL. 

## Unknown Alpha

The codes in this folder create the result in the section 4.2 in the paper. In this section, we don't know the true value of (alpha_0^-, alpha_0^+, alpha_1^-, alpha_1^+). We will perform on a sensitivity analysis on (alpha_1^-, alpha_1^+). For the value of alpha_0's, we assume two cases. The first case, we correctly specify the value of alpha_0. In the second case, we incorrectly specify the value of alpha_0. 

NOTE: the R file with "v2.R" is case where we misspecify the value of alpha_0.

These codes are written to run on Bluehive to create a heatmap. When we combine with the batch codes, we can create 500 heatmap with the sensitivity. 

The function to generate data set is included in each file. 

"AllCorrect.R" and "AllCorrectv2.R" correspond to case 1 in the overview section. Again, "...v2.R" file is where we misspecify the value of alpha_0. 
"FzMis.R" and "FzMisv2.R" correspond to case 2 in the overview section. Again, "...v2.R" file is where we misspecify the value of alpha_0.
"FazMis.R" and "FazMisv2.R" correspond to case 3 in the overview section. Again, "...v2.R" file is where we misspecify the value of alpha_0.
"Qmis.R" and "Qmisv2.R" correspond to case 3 in the overview section. Again, "...v2.R" file is where we misspecify the value of alpha_0.

## Heat Map Generator
The codes from from Unknow Alpha folders will give us 500 csv files for each case. The "Heat Map Generator all plots.R" will take those files and create all the heat maps in the section 4.2 in the paper. 

## Data Application

"Data Analysis.Rmd" import, clean the data and apply the proposed methods, OWL method, and Eric's method to this data set. It will create the heat maps in the section 5 in the paper. 

The data application can be found at the ENGAGE folder in Box.
