# modestum_divergence_wgs
This repository will contain details and data processing steps used in Amoroso et al. (in review). The details and data processing steps will be released upon acceptance of the manuscript. Please email Keaka Farleigh (keakafarleigh@virginia.edu; keakafarleigh@gmail.com) if you have any questions.

## Citation

The citation will be here.

## Genomic data

The genomic data for this project are deposited at [PRJNA1395130](); the link will not work until the BioProject is published.

## Structure of this repository.

This repository contains the computational workflow for [Amoroso et al. (in review)](). This ReadMe will tell you which scripts are associated with each analyses. Each script referenced herein is located in the `scripts/` directory.

### Data sequencing, assembly, and filtering

The original dataset was assembled and filtered by Novogene following [GATK's best practices workflow](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows). We then applied additional filters before inferring population structure and overall genetic differentiation. The additional filtering workflow is included in `ingroup_filtering.sh`. We generated a vcf file with invariant sites following the recommendations of pixy authors. This workflow can be found on the [pixy website](https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html).

### Population structure and demographic history

We evaluated the intraspecific structure of our dataset using the `admixture.sh` script and estimated phylogenies following the [captus tutorial](https://edgardomortiz.github.io/captus.docs/).

### Intraspecific variation

We estimated intraspecific differentiation between populations using Fst as estimated in [PopGenHelpR](https://kfarleigh.github.io/PopGenHelpR/). The script we used is named `differentiation.R`.

### Population divergence

We estimated population divergence and variation using pixy. This script is named `pixy.sh`. We identified and visualized candidate windows following the Irwin models using `candidate_windowidentification.R`. This script also contains code to perform Tajima's D simulations.

### Recombination rate


### Gene ontology
