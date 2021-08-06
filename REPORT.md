---
title: "REPORT"
author: "Bole and Oscar"
date: "8/2/2021"
output: .md document
---
# Genomic analysis reveals independent evolution of Plasmodium falciparum populations in Ethiopia; Reproducing methods
## Introduction 
_Plasmodium falciparum_ parasite populations in Ethiopia have been experiencing selective local pressures from drugs and immunity, leading to evolutionary adaptation. However, there was a paucity of data on genomic characterization and evolutionary adaptations of P. falciparum isolates from the central area of Ethiopia.
Whole-genome analysis of 25 P. falciparum isolates from central Ethiopia, specifically from West Arsi, were studied to determine their genetic diversity, population structures, and signatures of selection in known drug resistance alleles against global isolates from Cambodia, Thailand, DR Congo, and Malawi.
## Description of approach 
In a team of two, we reproduced the method of the above paper. Data used for reproducing was obtained from the link provided in the paper. The methods involved using GATK tools for variant discovery, snpEff tool for annotation, PLINK for population structure, and PCA.
The exercise was carried out using the HPC available at the centre. It is composed of a master node and two worker nodes: hpc01 and hpc02. Hpc01 has 256GB RAM, 22TB storage space, and 64 computer cores, while hpc02 contains 64GB RAM, 10TB storage, and 64 compute cores. For the exercise, the analysis was undertaken using hpc01. The data was stored in a shared accreditation folder, accessible to the members of the accreditation exercise.
We used Git and GitHub for collaboratively creating scripts, sharing results, and discussing the output. This tool allowed us to collaborate effectively.
## Obtaining data 
we obtained the data through the provided link in the paper[click here](https://www.malariagen.net/data/pf3k-5)
We downloaded the reference genome and already aligned Bam files. The script used to download the data. 
````
# downloading reference file
wget ftp://ngs.sanger.ac.uk/production/pf3k/release_5/Pfalciparum.genome.fasta.gz 
# downloading the bam files
wget ftp://ngs.sanger.ac.uk/production/pf3k/release_5/BAM/*.bam
````
## Indexing and creating a dictionary
Indexing a genome can be explained similar to indexing a book. If you want to know on which page a certain word appears or a chapter begins, it is much more efficient/faster to look it up in a pre-built index than going through every page of the book until you found it. Same goes for variant discovery. Indices allow, e.g. Haplotypecaller to narrow down the potential origin of a query variant within the file, saving both time and memory. Samtools index was used for this.
`CreateSequenceDictionary` creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by processing and analysis tools.
The script for indexing of the variants (BAM file) is;
````
#indexing of bam files
for file in *bam;
    do
    samtools index -b $file
done
````
and for creating dictionary and indexing of reference genome as below;
````

#creating dictionary for refrence genome
gatk CreateSequenceDictionary -R Pfalciparum.genome.fasta
#creating index for refrence genome
samtools faidx Pfalciparum.genome.fasta
````

## Variant calling
variant calling was done using GATK's Haplotypecaller "Best practices". Haplotypecaller calls germline SNPs and indels via local re-assembly of haplotypes. HaplotypeCaller runs per-sample to generate an intermediate GVCF, which can then be used in GenotypeGVCFs for joint genotyping of multiple samples in a very efficient way. In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. 
The script used for Haplotypecaller
```` 
#Hapolotype
parallel 'gatk HaplotypeCaller -R Pfalciparum.genome.fasta  -I {} -O {}.hppc.g.vcf' -ERC GVCF ::: *.bam

````  
and the result is ![here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Haplotypecaller.png).
## Combining.
We merged our HaplotypeCaller GVCF files into a single GVCF with appropriate annotations using CombineGVCFs.CombineGVCFs is meant to be used for merging of GVCFs that will eventually be input into GenotypeGVCFs.
The script used for CombineGVCFs
````
ls vcf* > vcf.list
gatk CombineGVCFs -R Pfalciparum.genome.fasta --variant vcfs.list -O combined.g.vcf.gz 
````
and result ![here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/combinedvcf.png).
## Genotyping.
We genotyped the combined GVCF file obtained from the CombineGVCFs tool using GenotypeGVCFs.This tool is designed to perform joint genotyping on a single input, containing one or many samples. In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with _-ERC GVCF_ or _-ERC BP_RESOLUTION_
The script used;
````
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R Pfalciparum.genome.fasta \
   -V combined.g.vcf.gz \
   -O genotyped.g.vcf.gz
````
and result ![here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Genotyped.png).
## Assigning a quality score
We used Hard-filtering germline short variants guidelines for filtering.
Hard-filtering consists of choosing specific thresholds for one or more annotations and throwing out any variants with annotation values above or below the set thresholds. By annotations, we mean properties or statistics describing each variant(QD, MQ, FS, DP, SOR,). 
### QD(QualByDepth)
This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. For filtering purposes, it is better to use QD than either QUAL or DP directly. There were two peaks where most variants are (around QD = 4 and QD = 28). These two peaks correspond to variants that are mostly observed in heterozygous (het) versus mostly homozygous-variant (hom-var) states, respectively, in the called samples. This is because hom-var samples contribute twice as many reads supporting the variant as do het variants.
### MQ(RMSMappingQuality)
This is the root mean square mapping quality over all the reads at the site. Instead of the average mapping quality of the site, this annotation gives the square root of the average of the squares of the mapping qualities at the site. It is meant to include the standard deviation of the mapping qualities. Including the standard deviation allows us to include the variation in the dataset. A low standard deviation means the values are close to the mean, whereas a high standard deviation means the values are far from the mean. When the mapping qualities are good at a site, the MQ will be around 60.

### FS(FisherStrand)
This is the Phred-scaled probability that there is strand bias at the site. Strand Bias tells us whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele. When there is little to no strand bias at the site, the FS value will be close to 0.

### SOR(StrandOddsRatio)
This is another way to estimate strand bias using a test similar to the symmetric odds ratio test. SOR was created because FS tends to penalize variants that occur at the ends of exons. Reads at the ends of exons tend to only be covered by reads in one direction, and FS gives those variants a bad score. SOR will take into account the ratios of reads that cover both alleles.

## Filtering
The filtering thresholds were obtained from the ggplot(data visualizing tool) on R-studio using a table generated from gatk's VariantsTotable (this tool extracts specified fields for each variant in a VCF file to a tab-delimited table)
The script of VariantsTotable
````gatk VariantsToTable -V genotyped.g.vcf.gz -F QD -F FS -F SOR -F MQ -F DP -O Hardfilter.table ```` 
and ggplots scripts are 
````
library(ggplot2)
QD.plot <- ggplot(data = Hardfilter, aes(x=QD)) + geom_density(alpha=0.2)
QD.plot
#generating MQ  PLOT
MQ.plot <- ggplot(data = Hardfilter, aes(x=MQ)) + geom_density(alpha=0.2)
MQ.plot
#Generating SOR plot
SOR.plot <- ggplot(data = Hardfilter, aes(x=SOR)) + geom_density(alpha=0.2)
SOR.plot
#Generating FS
FS.plot <- ggplot(data = Hardfilter, aes(x=FS)) + geom_density(alpha=0.2)
FS.plot + scale_x_log10()
````
and the results are [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Hardfilter.table) and result [here](https://github.com/bolekj/Plasmodium_falciparum/tree/master/ggplot/.
After obtaining the thresholds, gatk's VariantFiltration tool (which is designed for hard-filtering variant calls based on certain criteria) was used. The script used;
````
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
 ````
## selecting.
SelectVariants was used for selecting. This tool makes it possible to select a subset of variants based on various criteria in order to facilitate certain analyses. The criteria we used in this tool were for selecting variants that have passed all the filtering thresholds, and from the output of this, we selected SNPs. the scripts we used 
````
#selects the pass variants
 gatk SelectVariants \
     -R Pfalciparum.genome.fasta \
     -V filtered.g.vcf.gz \
     --select 'vc.isNotFiltered()' \
     -O selected.g.vcf.gz
````     
and
````
gatk SelectVariants \
    -V selected.g.vcf.gz \
    -select-type SNP \
    -O snps.g.vcf.gz
````    
The results on the same are ![here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Select.png) and ![here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Snps%20selected.png) respectively.
## snpEff
SnpEff is an open-source tool that annotates variants and predicts their effects on genes by using an interval forest approach. This program takes pre-determined variants listed in a data file that contains the nucleotide change and its position and predicts if the variants are deleterious.

### Building database.
SnpEff needs a database to perform genomic annotations. For our case it required us to build the database for Pfalciparum(3D7v.31) using snpEff.config (configures a new genome).We added genome entry to snpEff's configuration;

1.Get a GFF file 
````
mkdir path/to/snpEff/data/3D7v.31
cd path/to/snpEff/data/3D7v.31
wget ftp.sanger.ac.uk/pub/project/pathogens/gff3/2015-08/Pfalciparum.gff3.gz #(path provided for the data from the paper)
mv Pfalciparum.gff3.gz genes.gff.gz
````
2.Input reference sequence.
````
mkdir path/to/snpEff/data/genome
cp path/to/Pfalciparum.genome.fasta genome/
mv Pfalciparum.genome.fasta 3D7v.31.fasta
````
3.Add the new genome to the config file
````
nano snpEff.config
# Database for plasmodium falciparum 3D7v.31
         3D7v.31.genome : plasmodium falciparum 3D7v.31
         3D7v.31.reference : ftp://ngs.sanger.ac.uk/production/pf3k/release_5/Pfalciparum.genome.fasta.gz
````
4.Create database
````
cd /path/to/snpEff
java -jar snpEff.jar build -gff3 -v 3D7v.31
````
### Annotation.
Annotation is the process of identifying the locations of genes and all of the coding regions in a genome and determining what those genes do. The script we used is as below;
````
java -Xmx8g -jar snpEff.jar 3D7v.31 /opt/data/oscarmwaura/data/snps.g.vcf.gz  > snp.ann.g.vcf.gz
````
The snpEff annotation gives three files as output; vcf file ![click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Annotated_snps.png) , txt file[click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/snpEff_genes.txt) and html [click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/snpEff_summary.html)
## Population structure: PCA
We then investigated population structure using principal components analysis. Examing population structure can give us a great deal of insight into the history and origin of populations. To perform a PCA on our snp.ann.vcf data, we used plink -version (1.9). The following are the steps we took;
### 1.Linkage pruning
One of the major assumptions of PCA is that the data we use is independent, so as a first step, we need to prune our dataset of variants that are in linkage. We began by creating a directory where we carried out all the steps necessary for PCA analysis.
The script we used to do linkage pruning is
````
# perform linkage pruning - i.e. identify prune sites
plink --vcf snp.ann.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out cichlids
````
When complete, it will write out two files cichlids. prune.in [click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/plink_output/plink_pca_results/cichlids.prune.in) and cichlids.prune.out [click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/plink_output/plink_pca_results/cichlids.prune.out). The first one is a list of sites that fell below our linkage threshold - i.e. those we should retain. The other file is the opposite of this.
### 2.Perform a PCA
Next, we rerun plink with a few additional arguments to get it to conduct a PCA;
````
# create pca
plink --vcf snp.ann.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract cichlids.prune.in \
--make-bed --mind --pca --out cichlids_pca
````
The output gave us a series of new files which; three plink binary [bim](https://github.com/bolekj/Plasmodium_falciparum/blob/master/plink_output/plink_pca_results/cichlids_pca.bim),[bed](https://github.com/bolekj/Plasmodium_falciparum/blob/master/plink_output/plink_pca_results/cichlids_pca.bed),[fam](https://github.com/bolekj/Plasmodium_falciparum/blob/master/plink_output/plink_pca_results/cichlids_pca.fam) files and two PCA output cichlids.eigenval [click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/plink_output/plink_pca_results/cichlids_pca.eigenval) cichlids.eigenvec [click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/plink_output/plink_pca_results/cichlids_pca.eigenvec)
### 3.Plotting the PCA output
We then turned to R to plot the analysis we have produced!
We first moved the plink output (the two PCA outputs) into the working directory and then loaded the tidyverse package.

We then used a combination of readr and the standard scan function to read in the data.
#### Cleaning up the data
The output we got after reading in and scanning were not in reasonable form, so we had to clean it up by first removing a nuisance column then gave our pca data.frame proper column names. We did this using the R version of grep. We then use paste0 to combine the columns.

With these variables created, we remade our data.frame. Note the use of as.tibble to ensure that we make a tibble for easy summaries.

#### Plotting the data
We first made a plot of the eigenvalues[click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/plink_output/plink_pca_results/Rplot_pca_pva.pdf)
We then calculated the cumulative sum of the percentage variance.
Next, we moved on to the actual plotting of our PCA. 
The Rscript we used for performing all this is as provided 
````
# load tidyverse package
library(tidyverse)
# read in data
pca <- read_table("./cichlids_pca.eigenvec", col_names = FALSE)
eigenval <- scan("./cichlids_pca.eigenval")

#CLEANING UP THE DATA
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "samples"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
#sorting pops
# location
pop <- rep(NA, length(pca$samples))
pop[grep("PA", pca$samples)] <- "Cambodia"
pop[grep("PD0504-C", pca$samples)] <- "Thailand"
pop[grep("PD0519-C", pca$samples)] <- "Thailand"
pop[grep("GB4", pca$samples)] <- "Congo"

# remake data.frame
pca <- as_tibble(data.frame(pca, pop))
#Remove unwanted columns
pca <- pca[,-22]

# first convert to percentage variance explained
pve <- data.frame(eigenvals = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(eigenvals, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = pop )) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue" ,"green"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
```` 
The pca plot we got is [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/plink_output/plink_pca_results/Rplot_pca.pdf)
## Conclusion
We successfully managed to do variant discovery, annotation of the variants and population structure using principal component analysis.
The results we obtained in all the steps are not as obtained in the paper because we did not get all the raw data when downloading.
We continued with the steps in the paper because the main aim of this exercise was to learn how to go around when given a project(see the workflow).

