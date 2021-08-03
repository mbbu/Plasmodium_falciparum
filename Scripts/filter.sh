gatk VariantFiltration \
   -R Pfalciparum.genome.fasta \
   -V genotyped.g.vcf.gz \
   -O filtered.g.vcf.gz\
   --filter-name "QD37"\
   --filter-expression "QD > 37.0" \
   --filter-name "FS60"\
   --filter-expression "FS > 60.0" \
   --filter-name "SOR5"\
   --filter-expression "SOR > 5.0" \
   --filter-name "MQ25"\
   --filter-expression "MQ < 25.0"\

