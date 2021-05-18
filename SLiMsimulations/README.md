# SLiM and R code supporting individual-based simulations of "Allorecognition genes drive reproductive isolation in Podospora anserina"

**by Ivain Martinossi**


-------------------------------------------------------------------------------
## Context
Individual based simulation of two loci interacting in Podospora anserina are run using the software SLiM.
R is used to run repeated simulations with SLiM and extract and plot data. Both R and SLiM code are provided here.
In a first scenario, the Vr genotype invades a V1R population. In a second scenario, V1R invades a Vr population.
See article for more details.

-------------------------------------------------------------------------------
## Content
For scenario 1, Vr invades V1R:
Vr_invades_V1R.R
vr_invades_V1R.slim

For scenario 2, V1R invades Vr:
V1R_invades_Vr.slim
V1R_invades_Vr.slim

-------------------------------------------------------------------------------
## How to use
Install SLiM and R. 
Modify the R code to add the correct path for the SLiM software, and the correct path for the SLiM script on your own computer.
Modify the path in the SLiM code to output the data in the desired directory.
Run R code.

-------------------------------------------------------------------------------
## What does it do?
Almost every line of code in R and SLiM codes have comments describing their detailed functions.
To give a general idea, the SLiM codes runs 1500 generation of a population, starting at one genotype depending on the scenario.
At generation 500, the selfing rate is adjusted to the desired level.
At generation 600, the alternative genotype is added, representing 1% of the population, simulating a migration event.
The R code allows to run the above simulation across different parameter combinations and to automatically manage the resulting datasets.
Between generation 900 and 1000, genotype frequencies are sampled and used to produce plots with R.