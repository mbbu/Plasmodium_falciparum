 #selects the pass variants
 gatk SelectVariants \
     -R Pfalciparum.genome.fasta \
     -V filtered.g.vcf.gz \
     --select 'vc.isNotFiltered()' \
     -O selected.g.vcf.gz

