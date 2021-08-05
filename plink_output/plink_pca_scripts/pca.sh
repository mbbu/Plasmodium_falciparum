# prune and create pca
plink --vcf snp.ann.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract cichlids.prune.in \
--make-bed --mind --pca --out cichlids_pca
