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



