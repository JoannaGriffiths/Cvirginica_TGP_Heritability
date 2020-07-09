# Data and code from Griffiths et al. (in review) Transgenerational plasticity and the capacity to adapt to low salinity in the eastern oyster, Crassostrea virginica

## Size Analyses
*Script*: Oyster_size_ANOVA.R

*Input files*: oyster_size.txt, parent_oyster_size.txt

*Description*: Script contains code for normalizing size data for larvae and parental adults and performing an ANOVA and a Tukey’s post-hoc test


## Mortality Analyses
*Script*: Oyster_mortality_ANOVA.R

*Input files*: oyster_mortality.txt

*Description*: Script contains code for mortality data for each cross and performing an ANVOA

## Biochem Analyses
*Script*: biochem_ANOVA.R

*Input files*: biochem_data.txt, oyster_size_eggQual.txt

*Description*: Script contains code for ANOVAs to determine whether egg quality differs among dam acclimation site, whether egg quality is related to a dam's successful spawn, and whether egg quality is related to size of larvae.

## Heritability Analyses
*Script*: Oyster_MCMCglmm.R

*Input files*: oyster_ped.txt, oyster_size.txt

*Description*: Script contains different models for estimating heritability and variance components for larval size. Oyster_ped.txt contains pedigree information for each family.
