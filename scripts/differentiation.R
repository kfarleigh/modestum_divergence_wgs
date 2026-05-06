# Purpose: Estimate individual heterozygosity and genetic differentiation of Phrynosoma modestum
# Author: Keaka Farleigh
# Email: keakafarleigh@virginia.edu
# Date: January 13th, 2025

### Load packages
library(PopGenHelpR)
library(vcfR)
library(hierfstat)
library(adegenet)

##############################################
##### Read in population assignment file #####
##############################################

Pops <- read.csv('../Metadata/modestum_samples_mdata_PGH.csv')

# Remove RRM_2081
Pops <- Pops[which(Pops$Sample != "RRM_2081"),]

# Check to make sure that the population map is order correctly
vcf <- read.vcfR("../vcf/modestum_ingroup_0.9missing_noRRM2081.recode.vcf")

gt <- extract.gt(vcf)

inds <- colnames(gt)

Pops$Sample %in% inds

Pops$Sample == inds

# clear out our space
remove(vcf, gt, inds)

gc()

###########################
##### Differentiation #####
###########################

Dif <- Differentiation(data = '../vcf/modestum_ingroup_0.9missing_noRRM2081.recode.vcf', pops = Pops)
