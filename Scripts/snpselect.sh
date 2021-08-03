gatk SelectVariants \
    -V selected.g.vcf.gz \
    -select-type SNP \
    -O snps.g.vcf.gz

