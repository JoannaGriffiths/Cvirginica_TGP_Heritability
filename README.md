# Data and code from Griffiths et al. (in review) Transgenerational plasticity and the capacity to adapt to low salinity in the eastern oyster, Crassostrea virginica

## Size Analyses
*Script*: Oyster_size_ANOVA.R

*Input files*: oyster_size.txt, parent_oyster_size.txt

*Description*: Script contains code for normalizing size data for larvae and parental adults and performing an ANOVA and a Tukey’s post-hoc test


## Mortality Analyses
*Script*: Oyster_mortality_ANOVA.R

*Input files*: oyster_mortality.txt

*Description*: Script contains code for mortality data for each cross and performing an ANVOA

*Script*: parent_mort.R

*Input files*: parent_mort.txt, cumulative_parent_mort.txt

*Description*: Script contains code for mortality data for adults outplanted at the two field acclimation sites. parent_mort.txt contains mrotality data after two years of outplant and the cumulative_parent_mort.txt contains the mortality collected every two months for adults outplanted at the field sites.

## Biochem Analyses
*Script*: biochem_ANOVA.R

*Input files*: biochem_data.txt, oyster_size_eggQual.txt

*Description*: Script contains code for ANOVAs to determine whether egg quality differs among dam acclimation site, whether egg quality is related to a dam's successful spawn, and whether egg quality is related to size of larvae.

## Heritability Analyses
*Script*: Oyster_MCMCglmm.R

*Input files*: oyster_ped.txt, oyster_size.txt

*Description*: Script contains different models for estimating heritability and variance components for larval size. Oyster_ped.txt contains pedigree information for each family. This script also uses the breeder's equation to estimate the required mortality to return larval sizes back to ambient salinity conditions after one generation of selection.

## Genetic Co-Variation Analyses
*script*: Genetic_coVariation.R

*Input files*: family_means_v2.csv, oyster_ped_covar_matrix.txt

*Description*: Script contains model for estimating genetic co-variation between families reared at low and high salinity conditions. oyster_ped_covar_matrix.txt contains pedigree information for each family. family_means_v2.csv contains mean larval size for each family reared at low and high salinity. 
