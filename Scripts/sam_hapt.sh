#indexing of bam files
for file in *bam;
    do
    samtools index -b $file
done

#Hapolotype
parallel 'gatk HaplotypeCaller -R Pfalciparum.genome.fasta  -I {} -O {}.hppc.g.vcf' -ERC GVCF ::: *.bam


