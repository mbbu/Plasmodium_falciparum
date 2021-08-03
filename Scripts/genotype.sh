 gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R Pfalciparum.genome.fasta \
   -V combined.g.vcf.gz \
   -O genotyped.g.vcf.gz


