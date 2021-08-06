---
title: "REPORT"
author: "Bole and Oscar"
date: "8/2/2021"
output: .md document
---
# Genomic analysis reveals independent evolution of Plasmodium falciparum populations in Ethiopia; Reproducing methods
## Introduction 
_Plasmodium falciparum_ parasite populations in Ethiopia have been experiencing local selective pressures from drugs and immunity, leading to evolutionary adaptation. However, there was a paucity of data on genomic characterization and evolutionary adaptations of P. falciparum isolates from the central area of Ethiopia.
Whole-genome analysis of 25 P. falciparum isolates from central Ethiopia, specifically from West Arsi, were studied to determine their genetic diversity, population structures, and signatures of selection in known drug resistance alleles against global isolates from Cambodia, Thailand, DR Congo, and Malawi.
## Description of approach 
In a team of two,we reproduced the method of the above paper.Data used for reproducing was obtained from the link provide in the paper.The methods involved using of gatk tools for variant discovery and snpEff tool for annotation.
The exercise was carried out using the HPC available at the center. It is composed of a master node and two worker nodes: hpc01 and hpc02. Hpc01 has 256GB RAM, 22TB storage space, and 64 computer cores, while hpc02 contains 64GB RAM, 10TB storage, and 64 compute cores. For the exercise, the analysis was undertaken using hpc01. The data was stored in a shared accreditation folder, accessible to the members of the accreditation exercise.
We used Git and GitHub for collaboratively creating scripts, sharing results, and discussing the output. This tool allowed us to collaborate effectively.
## Obtaining data 
we obtained the data through the provided link in the paper[click here](https://www.malariagen.net/data/pf3k-5)
We downloaded reference genome and already aligned Bam fileS. The script used to download the data is ```` wget ftp://ngs.sanger.ac.uk/production/pf3k/release_5/Pfalciparum.genome.fasta.gz ```` and ````wget ftp://ngs.sanger.ac.uk/production/pf3k/release_5/BAM/*.bam```` respectively.
## Indexing and creating dictionary
Indexing a genome can be explained similar to indexing a book. If you want to know on which page a certain word appears or a chapter begins, it is much more efficient/faster to look it up in a pre-built index than going through every page of the book until you found it. Same goes for variant discovery. Indices allows,e.g Haplotypecaller  to narrow down the potential origin of a query variant within the file, saving both time and memory.
Creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by processing and analysis tools.
The script for indexing of the variants (BAM file) is [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Scripts/sam_hapt.sh) and for creating dictionary and indexing of refrence genome as below;
````
#creating dictionary for refrence genome
gatk CreateSequenceDictionary -R Pfalciparum.genome.fasta
#creating index for refrence genome
samtools faidx Pfalciparum.genome.fasta
````

## Variant calling
varint calling was done using gatk's Haplotypecaller "Best practices". Haplotypecaller calls germline SNPs and indels via local re-assembly of haplotypes. HaplotypeCaller runs per-sample to generate an intermediate GVCF, which can then be used in GenotypeGVCFs for joint genotyping of multiple samples in a very efficient way. In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. 
The script used for Haplotypecaller and result are [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Scripts/sam_hapt.sh) and [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Haplotypecaller.png) respectively.
## Combining.
We merged our HaplotypeCaller GVCF files into a single GVCF with appropriate annotations using CombineGVCFs.CombineGVCFs is meant to be used for merging of GVCFs that will eventually be input into GenotypeGVCFs.
The script used for CombineGVCFs and result are [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Scripts/combined.sh) and [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/combinedvcf.png) respectively.
## Genotyping.
We genotyped the combined GVCF file obtained from the CombineGVCFs tool using GenotypeGVCFs.This tool is designed to perform joint genotyping on a single input, which may contain one or many samples. In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with _-ERC GVCF_ or _-ERC BP_RESOLUTION_
The script used for GenotypeGVCFs and result are [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Scripts/genotype.sh) and [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Genotyped.png) respectively.
## Assigning quality score
We used Hard-filtering germline short variants guidelines for filtering.
Hard-filtering consists of choosing specific thresholds for one or more annotations and throwing out any variants that have annotation values above or below the set thresholds. By annotations, we mean properties or statistics that describe for each variant(QD,MQ,FS,DP,SOR,). 
### QD(QualByDepth)
This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. For filtering purposes it is better to use QD than either QUAL or DP directly. There was two peaks where the majority of variants are (around QD = 4 and QD = 28). These two peaks correspond to variants that are mostly observed in heterozygous (het) versus mostly homozygous-variant (hom-var) states, respectively, in the called samples. This is because hom-var samples contribute twice as many reads supporting the variant than do het variants.
### MQ(RMSMappingQuality)
This is the root mean square mapping quality over all the reads at the site. Instead of the average mapping quality of the site, this annotation gives the square root of the average of the squares of the mapping qualities at the site. It is meant to include the standard deviation of the mapping qualities. Including the standard deviation allows us to include the variation in the dataset. A low standard deviation means the values are all close to the mean, whereas a high standard deviation means the values are all far from the mean.When the mapping qualities are good at a site, the MQ will be around 60.
### FS(FisherStrand)
This is the Phred-scaled probability that there is strand bias at the site. Strand Bias tells us whether the alternate allele was seen more or less often on the forward or reverse strand than the reference allele. When there little to no strand bias at the site, the FS value will be close to 0.
### SOR(StrandOddsRatio)
This is another way to estimate strand bias using a test similar to the symmetric odds ratio test. SOR was created because FS tends to penalize variants that occur at the ends of exons. Reads at the ends of exons tend to only be covered by reads in one direction and FS gives those variants a bad score. SOR will take into account the ratios of reads that cover both alleles.
## Filtering
The filtering thresholds were obtained from the ggplot(data visualizing tool) on R-studio using a table generated from gatk's VariantsTotable (this tool extracts specified fields for each variant in a VCF file to a tab-delimited table)
The script of VariantsTotable and ggplots are [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Scripts/table.sh) and [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Scripts/ggplot_scripts.R) and the results are [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Hardfilter.table) and [here](https://github.com/bolekj/Plasmodium_falciparum/tree/master/ggplots) respectively.
After obtaining the thresholds, gatk's VariantFiltration tool (which is designed for hard-filtering variant calls based on certain criteria) was used.The script used is [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Scripts/filter.sh) 
## selecting.
SelectVariants was used for selecting. This tool makes it possible to select a subset of variants based on various criteria in order to facilitate certain analyses. The criteria we used in this tool were for selecting variants that have passed all the filtering threshholds and from the output of this we selected SNPs. the scripts we used are [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Scripts/select.sh) and [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Scripts/snpselect.sh) respectively.The results on the same are [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Select.png) and [here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Snps%20selected.png) respectively.
## snpEff
SnpEff is an open source tool that annotates variants and predicts their effects on genes by using an interval forest approach. This program takes pre-determined variants listed in a data file that contains the nucleotide change and its position and predicts if the variants are deleterious.
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
Annotation is the process of identifying the locations of genes and all of the coding regions in a genome and determining what those genes do.The script we used is as below;
````
java -Xmx8g -jar snpEff.jar 3D7v.31 /opt/data/oscarmwaura/data/snps.g.vcf.gz  > snp.ann.g.vcf.gz
````
The snpEff annotation gives three files as output; vcf file[click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/Results_snapshots/Annotated_snps.png) , txt file[click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/snpEff_genes.txt) and html [click here](https://github.com/bolekj/Plasmodium_falciparum/blob/master/snpEff_summary.html)
## Population structure: PCA
We then investigated population structure using principal components analysis.Examing population structure can give us a great deal of insight into the history and origin of populations. To perform a PCA on our cichlids data, we used plink -version (1.9). The following are the steps we took;
### 1.Linkage pruning
One of the major assumptions of PCA is that the data we use is indpendent, so as a first step, we need to prune our dataset of variants that are in linkage. We begun by creating a directory were we carried out all the steps necessary for pca analysis.
The script we used to do linkage pruning is [here]().When complete, it will write out two files cichlids.prune.in[click here]() and cichlids.prune.out[click here](). The first one is a list of sites which fell below our linkage threshold - i.e. those we should retain. The other file is the opposite of this.
### 2.Perform a PCA
Next we rerun plink with a few additional arguments to get it to conduct a PCA.[click here]() for the script we used.The output gave us a series of new files which; three plink binary [bim](),[bed](),[fam]() files and two PCA output cichlids.eigenval[click here]() cichlids.eigenvec[click here]()
### 3.Plotting the PCA output
We then turned to R to plot the analysis we have produced!
We first moved the plink output (the two PCA output) into the working directory,then loaded the tidyverse package.
We then used a combination of readr and the standard scan function to read in the data.
#### Cleaning up the data
The output we got after reading in and scanning  was not in reasonable form, so we had to clean it up by first removing a nuisance column then gave our pca data.frame proper column names.
. We will do this using the R version of grep. We then use paste0 to combine the columns.
With these variables created, we remade our data.frame. Note the use of as.tibble to ensure that we make a tibble for easy summaries
#### Plotting the data
We first made a plot of the eigenvalues. It is quite straightforward to translate these into percentage variance explained. We then created bar plot showing the percentage of variance each principal component explains.
We could calculate this with the cumsum function
Next we moved on to the actual plotting of our PCA. 
The Rscript we used for performing all this is as provided [here]() 
## Conclusion
We successfully managed to do variant discovery and annotation of the variants.

