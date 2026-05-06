# The all sites vcf was generated follwing https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html
bcftools index -s allsites_final.vcf.gz | cut -f 1 | while read C; do bcftools view -O z -o split.${C}.vcf.gz allsites_final.vcf.gz "${C}" ; done

mkdir ./Pixy_output

conda activate pixypy #python environment with pixy, htslib, samtools installed

#run pixy, pixy_populations is a text file containing population assignment for each individual (K=2)
for i in `ls split*`; do tabix $i; pixy --stats pi fst dxy --vcf $i --populations pixy_populations.txt --window_size 10000 --n_cores 8 --output_folder ./Pixy_output --output_prefix $i; done
