# Purpose: This script runs ADMIXTURE
# Author: Julia Amoroso
# Date: November 3rd, 2024
# Email: jamoroso@gradcenter.cuny.edu

#generate ADMIXTURE input file
plink --vcf rn_modestum_ingroup_thinned.vcf.gz --make-bed --out admixture_modestum --allow-extra-chr

#rename chromosomes to ADMIXTURE format
awk '{$1=0; print $0}' admixture_modestum.bim > admixture_modestum.bim.tmp
mv admixture_modestum.bim.tmp admixture_modestum.bim

#run admixture for K = 1 through K = 10
for i in {1..10}
do
  admixture --cv admixture_modestum.bed $i > log${i}.out
done
