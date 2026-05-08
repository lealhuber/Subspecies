### PCA analysis
## adapted from https://speciationgenomics.github.io/pca/

library(tidyverse)

datdir <- "/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/Subspecies/data/PCA_prelim/"

## read in data
pca <- read_table(paste0(datdir,"pca_all.eigenvec"), col_names = FALSE)
eigenval <- scan(paste0(datdir,"pca_all.eigenval"))

## clean up
pca <- pca[,-1] # remove double id column
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
## add species and sex info
specvec <- c(rep("black",10), rep("blue",8), rep("red",8))
pca$pop <- specvec
pca$sex <- c(rep("M",5),rep("F",5),rep("M",3), rep("F",5),rep("M",4),rep("F",4))


### plotting
# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# barplot of how much explained variance per pc
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
# there will actually be an argument to just plot the first 2
cumsum(pve$pve) # calculate the cumulative sum of the percentage variance explained
# first 2 explain almost a quarter
## plot pca
library(ggrepel)
ggplot(pca, aes(PC1, PC2, col = pop, shape = sex, label = ind)) + geom_point(size = 3) +
  scale_colour_manual(values = c("black", "blue", "red"))+
  geom_text_repel()+
  coord_equal() + theme_light(base_size = 17)+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
# PC1 separates red from black and blue, PC2 separates black from blue
# 133 and 122 are both between black and blue, so they are probably both half each
# Thus more evidence that 133 is not red at all

## Plot the PC histogram
barplot(pve$pve, names.arg = pve$PC, xlab = "Principal Component", ylab = "% of variance explained")

#### again with new data ---------------------------------------------------

datdir <- "/Users/Lea/Library/CloudStorage/OneDrive-Aarhusuniversitet/LeaH/Subspecies/data/PCA/"

## read in data
pca <- read_table(paste0(datdir,"pca_all.HWE.eigenvec"), col_names = FALSE)
eigenval <- scan(paste0(datdir,"pca_all.HWE.eigenval"))

## clean up
pca <- pca[,-1] # remove double id column
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
## add species and sex info
specvec <- c(rep("red",6), rep("black",10), rep("blue",8))
pca$pop <- specvec
pca$sex <- c(rep("M",4),"F","F", rep("M",5), rep("F",5),rep("M",4),rep("F",4))


### plotting
# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
# barplot of how much explained variance per pc
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
# there will actually be an argument to just plot the first 2
cumsum(pve$pve) # calculate the cumulative sum of the percentage variance explained
# first 2 explain almost a quarter
## plot pca
library(ggrepel)
ggplot(pca, aes(PC1, PC2, col = pop, shape = sex, label = ind)) + geom_point(size = 3) +
  scale_colour_manual(values = c("black", "blue", "red"))+
  geom_text_repel()+
  coord_equal() + theme_light(base_size = 17)+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## Plot the PC histogram
barplot(pve$pve, names.arg = pve$PC, xlab = "Principal Component", ylab = "% of variance explained (eigenvalues)")



