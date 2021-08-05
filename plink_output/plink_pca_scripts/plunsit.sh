 # perform linkage pruning - i.e. identify prune sites
plink --vcf snp.ann.vcf --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out cichlids
