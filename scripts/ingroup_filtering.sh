# Purpose: This script prepares an ingroup vcf for ADMIXTURE and PopGenHelpR
# Author: Julia Amoroso
# Date: November 3rd, 2024
# Email: jamoroso@gradcenter.cuny.edu


#remove outgroup from vcf
vcftools --gzvcf merge.filter.snp.vcf.gz --remove-indv UWBM_7227 --recode --stdout | gzip -c > \modestum_ingroup.vcf.gz

#filter vcf for quality, missingness, linkage disequilibrium, make SNPs bi-allelic
vcftools --gzvcf modestum_ingroup.vcf.gz --minQ 20 --max-missing 0.9 --min-alleles 2 --max-alleles 2 --maf 0.05 --thin 10000 --recode --stdout | gzip -c > \modestum_ingroup_thinned.vcf.gz
