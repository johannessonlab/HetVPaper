#!/usr/bin/env Rscript

### SNPpop: Variant calling and population structure in *Podospora anserina* using SNPRelate
#############################################################################
# Produce PoCA of mating data and PCAs of SNPs from the Wageningen collection using SNPRelate
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-06-07, 2019-11-07
# Version 2
# =======================================
# http://corearray.sourceforge.net/tutorials/SNPRelate/
# ## Install the package from Bioconductor repository:
# source("http://bioconductor.org/biocLite.R")
# biocLite("gdsfmt")
# biocLite("SNPRelate")

# ============================
# Load the necessary libraries
# ============================
# Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(dplyr)
library(cluster)
library(vegan) #eigenvals
library(cowplot) # For plot_grid

# ============================
# Reading the data
# ============================
# Metadata
PopData <- read.csv(snakemake@input[[1]], header = TRUE, sep=",")
# PopData <- read.csv("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/data/SNPpopMetadata.csv", header=TRUE, sep=",")

# Load distance matrix
podandm <- read.csv(snakemake@input[[2]], na.strings = "na", sep = ",", header = TRUE)
# podandm <- read.csv(file="/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/data/PodosporaMating.csv", na.strings = "na", sep = ",", header = TRUE)

# SNP data
vcf.fn <- snakemake@input[[3]] # The vcf file
# vcf.fn <- "/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/results/PodoPop-snps-NoTEs-gatkPASS-NoSibs-miss1.vcf.gz"

workinggds <- snakemake@output[[1]] # For the produced gds file
# workinggds <- "/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/temp/PodoPop-snps-NoTEs-gatkPASS-NoSibs-miss1.gds"

cores = snakemake@threads
# cores = 4

# ============================
# PCoA of mating data
# ============================
# Remove little s because we didn't sequence that one
podandm <- podandm[podandm$X != "s", colnames(podandm) != "s"]

print("Calculating PCoA of the mating data from van der Gaag (2005) ...", quote=FALSE)

# The data is not actually numerical but ordinal
podandm[sapply(podandm, is.numeric)] <- lapply(podandm[sapply(podandm, is.numeric)], factor, order = TRUE, levels = c("1", "2", "3", "4", "5", "6")) 

# Convert distance matrix to dissimilarity matrix with gower distance (data is ordinal)
d2 <- daisy(podandm %>% select(-X), metric = "gower") # The column of strain ID has to be removed, otherwise it turns them into binary data and artificially adds them to the distances

# Use Partitioning Around Medioids (pam) to determine membership to K clusters
pam_fit2 <- pam(d2, k=2)

# ## PCA Plot of pam result
# clusplot(pam_fit2, labels=4, lines=2, color=TRUE, col.clus = c("coral3", "cadetblue3")) # It calls the object stored in the $call slot.
# # labels= 4, only the ellipses are labelled in the plot.
# # lines,	integer out of 0, 1, 2, used to obtain an idea of the distances between ellipses. 
# # lines = 2, a line segment between the boundaries of the distance between two ellipses E1 and E2 is drawn (along the line connecting the centers m1 and m2).

# # But this plot is very ugly, so I'll plot it using the clustering groups (pam_fit2$clustering)
# podandm %>% colnames()
# pam_fit2$clustering

## --- Gap statistic ---
set.seed(123)
gap_stat <- clusGap(d2 %>% as.matrix(), FUN = pam, nstart = 25,
                    K.max = 10, B = 500)
plot(gap_stat, main = "clusGap(., FUN = pam, nstart=25, B= 500)")

### ----- Recalculate the PCoA -----
### Getting ammount of variance explained
## I got this from https://github.com/wilkox/doPCoA/blob/master/R/do_PCoA.R
getvariancePCoA <- function(PCoA){
  # Calculate variance explained
  # Method is from http://r-sig-ecology.471788.n2.nabble.com/\
  # Variability-explanations-for-the-PCO-axes-as-in-Anderson-and-Willis-2003-td6429547.html
  Eigenvalues <- eigenvals(PCoA)
  Variance <- Eigenvalues / sum(Eigenvalues) 
  Variance1 <- (100 * signif(Variance[1], 4)) %>% round(digits=2)
  Variance2 <- (100 * signif(Variance[2], 4)) %>% round(digits=2)
  
  # ‘Variance1’ and ‘Variance2’ are the percentage variance explained by the first and second PCoA axes respectively.
  ## From the original source
  # Things can be a bit more complex with PCO as some dissimilarity 
  # coefficients are not metric and thus cannot be represented in an 
  # Euclidean space (what PCO is trying to do) so we get negative 
  # eigenvalues. The software can deal with this by adding a value to the 
  # dissimilarities to avoid the negative eigenvalues, but just be aware of 
  # this issue - the negative eigenvalues, if they occur, represent variance 
  # in imaginary dimensions (as we have imaginary numbers and normal numbers 
  # [I forget the correct term, rational numbers?] in mathematics) and 
  # possibly need to be considered. 
  
  return(c(Variance1, Variance2))
}

### Getting ammount of variance explained
# following http://r-sig-ecology.471788.n2.nabble.com/Variability-explanations-for-the-PCO-axes-as-in-Anderson-and-Willis-2003-td6429547.html 
# Because of the negative values of the eigenvalues
# cmdscale(d2, k=NROW(podandm %>% select(-X))-1, eig = TRUE)$eig  #only 31 of the first 45 eigenvalues are > 0
PCoA <- cmdscale(d2, k=NROW(podandm %>% select(-X))-1, eig = TRUE, add = TRUE)
varianza <- getvariancePCoA(cmdscale(d2, k=NROW(podandm %>% select(-X))-1, eig = TRUE, add = TRUE)) # 16.00  5.88 # This matches the output of clusplot perfectly!
# But the coordiantes are the same as cmdscale(d2, eig = TRUE)!

## Make a dataframe with all info for plotting
# Add the eigenvectors
matingpcoa <- data.frame(PCoA$points) #eig gives the eigenvalues
names(matingpcoa) <- c("PC1", "PC2")
mgsamples <- colnames(podandm)[-1]

matingpcoa <- cbind(SimpleID = mgsamples, Cluster = factor(pam_fit2$clustering), matingpcoa[c("PC1", "PC2")])
matingpcoa <- merge(matingpcoa, PopData[c("Sample", "SimpleID")], by = "SimpleID", all.x = TRUE) # Put the Samples genome names too

# # Index of little s
# little_s <- which(matingpcoa$SimpleID %in% "s")

# I don't need the SimpleID anymore
matingpcoa <- within(matingpcoa, rm("SimpleID")) # or snpschr5pca <- subset(snpschr5pca, select = -c(SimpleID))
# Add het-v just to check
matingpcoa <- merge(matingpcoa, PopData[c("Sample", "Het.v")], by = "Sample", all.x = TRUE)

# Make a theme for the size of labels
axessizetheme <- theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 16))

# The actual plot
matingplot <- ggplot(matingpcoa, aes(x = PC1, y = PC2, colour = Cluster, fill = Cluster)) +  geom_point(size = 5, alpha = 0.8) +
  stat_ellipse(alpha = 0.05, level = 0.99, geom = "polygon") +
  xlab(paste0("PC1 (", varianza[1], "%)")) + ylab(paste0("PC2 (", varianza[2], "%)")) +
  theme(panel.background = element_rect(fill = NA, color = "black"), # The fill = NA is necessary of it will plot over the points a white background
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20)) + # Title in the center 
  axessizetheme +
  ggtitle("Mating data") + scale_color_manual(values=c("#B66DE9", "darkseagreen3"))

# ============================
# Process SNPs
# ============================
# Reformat
snpgdsVCF2GDS(vcf.fn, workinggds, method="biallelic.only", ignore.chr.prefix="chromosome_") # By default the prefix is "Chr"

# Open the GDS file
genofile <- snpgdsOpen(workinggds)

# ============================
# PCA of whole genome SNPs
# ============================
print("Plotting a PCA of all the SNPs in the genome ...", quote=FALSE)

## LD-based SNP pruning
set.seed(1000)
# Try different LD thresholds for sensitivity analysis (and MAF)
snpset <- snpgdsLDpruning(genofile, ld.threshold=1, autosome.only=FALSE, maf = 0.01, num.thread = cores, slide.max.bp=10000) # Without the MAF filtering, S behaves weird; the value of slide.max.bp is 10 kb to match the genome scans and because LD decay is already very low at that point for most chromosomes. 

# Get all selected snp id
snpset.id <- unlist(snpset)

## PCA
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=FALSE, snp.id=snpset.id)

# variance proportion (%)
pc.percent <- pca$varprop*100

# Add the PCs into a df
allsnppca <- data.frame(Sample = pca$sample.id,
                        EV1 = pca$eigenvect[,1],    # the first eigenvector
                        EV2 = pca$eigenvect[,2],    # the second eigenvector
                        EV3 = pca$eigenvect[,3],    # the third eigenvector
                        stringsAsFactors = FALSE)
allsnppca <- merge(allsnppca, PopData, by = "Sample")

## Plot the mating group on top
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
# prepare the colors
clustV <- matingpcoa %>% filter(Cluster == 1) # Het-V
clustA <- matingpcoa %>% filter(Cluster == 2) # Het-V1

colors <- rep("gray", length(samples))
colors[samples %in% clustV$Sample]<- "#B66DE9"
colors[samples %in% clustA$Sample]<- "darkseagreen3"

# Add it to the dataframe
allsnppca <- cbind(allsnppca, colors)

# PC1 vs PC2
allsnppca12 <- ggplot(allsnppca, aes(x = EV1, y = EV2, colour = colors, fill = colors)) +  geom_point(size = 5, alpha = 0.8) +
  xlab(paste0("PC1 (", round(pc.percent[1], digits =2), "%)")) + ylab(paste0("PC2 (", round(pc.percent[2], digits =2), "%)")) +
  theme(panel.background = element_rect(fill = NA, color = "black"), # The fill = NA is necessary of it will plot over the points a white background
        legend.position = "none") + 
  scale_color_manual(values=c("#B66DE9", "darkseagreen3", "#999999"))

# PC1 vs PC3
allsnppca13 <- ggplot(allsnppca, aes(x = EV1, y = EV3, colour = colors, fill = colors)) +  geom_point(size = 5, alpha = 0.8) +
  xlab(paste0("PC1 (", round(pc.percent[1], digits =2), "%)")) + ylab(paste0("PC3 (", round(pc.percent[3], digits =2), "%)")) +
  theme(panel.background = element_rect(fill = NA, color = "black"), # The fill = NA is necessary of it will plot over the points a white background
        legend.position = "none") + 
  scale_color_manual(values=c("#B66DE9", "darkseagreen3", "#999999"))

# Plot the whole genome data
plot_row <- plot_grid(allsnppca12, allsnppca13, ncol=2)
title <- ggdraw() + draw_label("Whole genome SNP data", size = 20)
PCAallgenome <- plot_grid(title, plot_row,  ncol = 1, rel_heights = c(0.1, 1) )

ggsave(snakemake@output[[2]], plot=PCAallgenome, width = 7, height = 4.5)
# ggsave("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/results/figures/FigS4_PaPCA_all.pdf", plot=PCAallgenome, width = 5, height = 10)


# ============================
# Calculate the SNP correlations between eigenvactors and SNP genotypes
# ============================
print("Calculating the SNP correlations between eigenvactors and SNP genotypes ...", quote=FALSE)

## To calculate the SNP correlations between eigenvactors and SNP genotypes:
# Get chromosome index
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4, num.thread = cores)

savepar <- par(mfrow=c(2,1), mai=c(0.45, 0.55, 0.1, 0.25))
for (i in 1:2) {
  plot(abs(CORR$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste("PC", i), col=chr, pch="+")
}

# ggplot2 version
chr <- paste("Chr", chr, sep="")  # The package automatically removes "Chr" from the chromosome name
CORRdf <- data.frame(pos = pos, corr1 = abs(CORR$snpcorr[1,]), corr2 = abs(CORR$snpcorr[2,]), corr3 = abs(CORR$snpcorr[3,]), corr4 = abs(CORR$snpcorr[4,]), chr = chr)


corrpc1 <- ggplot(CORRdf, aes(x = pos, y = corr1, colour = chr)) + geom_point(size = 1, alpha = 0.8) +
  facet_grid(chr ~ .) + ylim(0, 1) + xlab("Chromosomal position (bp)") + ylab("Correlation coefficient") +
  ggtitle("PC1") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA, colour = "grey20"), 
        panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85", colour = "grey20"))

corrpc2 <- ggplot(CORRdf, aes(x = pos, y = corr2, colour = chr)) + geom_point(size = 1) +
  facet_grid(chr ~ .) + ylim(0, 1) + xlab("Chromosomal position (bp)") + ylab("Correlation coefficient") +
  ggtitle("PC2") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.background = element_rect(fill = "white", colour = NA), 
        panel.border = element_rect(fill = NA, colour = "grey20"), 
        panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85", colour = "grey20"))

pdf(snakemake@output[[3]], height=7, width=15)
# pdf("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/results/Figures/PaPCA_corr.pdf", height=7, width=13)
plot_grid(corrpc1, corrpc2, ncol=2)
dev.off()

# ============================
# PCA of SNPs by chromosome
# ============================
print("Plotting a SNPs PCA per chromosome ...", quote=FALSE)

snpid <- read.gdsn(index.gdsn(genofile, "snp.id"))
CORRdf <- cbind(snpid, CORRdf)

## Function to plot the PCA of each chromosome for different combinations of eigenvalues
pcaperchr <- function(genofile, chr = 1, cores = 4, LD = 1, pc1 = 1, pc2 = 2){
  chrname = paste0("Chr", chr)
  
  # Clean small frequency variants
  snpset <- snpgdsLDpruning(genofile, ld.threshold=LD, autosome.only=FALSE, maf = 0.01, snp.id= CORRdf %>% filter(chr == chrname) %>% .$snpid)
  # Get all selected snp id
  snpset.id <- unlist(snpset)
  # Do PCA with subset
  pcachr <- snpgdsPCA(genofile, num.thread=cores, autosome.only=FALSE, snp.id= snpset.id)
  
  # variance proportion (%)
  pc.percent <- pcachr$varprop*100
  
  # Add the PCs into a df
  tab <- data.frame(Sample = pcachr$sample.id,
                    EV1 = pcachr$eigenvect[,pc1],    # the first eigenvector
                    EV2 = pcachr$eigenvect[,pc2],    # the second eigenvector
                    stringsAsFactors = FALSE)
  tab <- merge(tab, PopData, by = "Sample")
  
  ## Plot the mating group on top
  samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
  # prepare the colors
  clustV <- matingpcoa %>% filter(Cluster == 1) # Het-V
  clustA <- matingpcoa %>% filter(Cluster == 2) # Het-V1
  
  colors <- rep("gray", length(samples))
  colors[samples %in% clustV$Sample]<- "#B66DE9"
  colors[samples %in% clustA$Sample]<- "darkseagreen3"
  
  # Add it to the dataframe
  tab <- cbind(tab, colors)
  
  # Mating group 1 is Het-V, and mating group 2 is Het-V1
  p <- ggplot(tab, aes(x = EV1, y = EV2, colour = colors, fill = colors)) +  geom_point(size = 5, alpha = 0.8) +
    xlab(paste0("PC", pc1, " (", round(pc.percent[pc1], digits=2), "%)")) + ylab(paste0("PC", pc2, " (", round(pc.percent[pc2], digits=2), "%)")) +
    theme(panel.background = element_rect(fill = NA, color = "black"), # The fill = NA is necessary of it will plot over the points a white background
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none") + 
    ggtitle(paste0("Chromosome ", chr)) + scale_color_manual(values=c("#B66DE9", "darkseagreen3", "#999999"))
  return(p)
}


PaPCA12 <- plot_grid(pcaperchr(genofile, chr = "1"),
                     pcaperchr(genofile, chr = "2"),
                     pcaperchr(genofile, chr = "3"),
                     pcaperchr(genofile, chr = "4"),
                     pcaperchr(genofile, chr = "5"),
                     pcaperchr(genofile, chr = "6"),
                     pcaperchr(genofile, chr = "7"),
                     ncol=4)
ggsave(snakemake@output[[4]], plot=PaPCA12, width = 13, height = 7)
# pdf("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/results/figures/PaPCA12.pdf", width = 15, height = 8)

PaPCA13 <- plot_grid(pcaperchr(genofile, chr = "1", pc2 = 3),
                     pcaperchr(genofile, chr = "2", pc2 = 3),
                     pcaperchr(genofile, chr = "3", pc2 = 3),
                     pcaperchr(genofile, chr = "4", pc2 = 3),
                     pcaperchr(genofile, chr = "5", pc2 = 3),
                     pcaperchr(genofile, chr = "6", pc2 = 3),
                     pcaperchr(genofile, chr = "7", pc2 = 3),
                     ncol=4)
ggsave(snakemake@output[[5]], plot=PaPCA13, width = 13, height = 7)
# pdf("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/results/figures/PaPCA13.pdf", width = 15, height = 8)

PaPCA23 <- plot_grid(pcaperchr(genofile, chr = "1", pc1 = 2, pc2 = 3),
                     pcaperchr(genofile, chr = "2", pc1 = 2, pc2 = 3),
                     pcaperchr(genofile, chr = "3", pc1 = 2, pc2 = 3),
                     pcaperchr(genofile, chr = "4", pc1 = 2, pc2 = 3),
                     pcaperchr(genofile, chr = "5", pc1 = 2, pc2 = 3),
                     pcaperchr(genofile, chr = "6", pc1 = 2, pc2 = 3),
                     pcaperchr(genofile, chr = "7", pc1 = 2, pc2 = 3),
                     ncol=4)
ggsave(snakemake@output[[6]], plot=PaPCA23, width = 13, height = 7)
# pdf("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/results/figures/PaPCA23.pdf", width = 15, height = 8)

# ============================
# PCoA of mating data and PCA of chr5
# ============================
chr5withline <- pcaperchr(genofile, chr = "5") + 
  geom_vline(xintercept = 0, col = "gray") +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=16)) # Adjust the size of the axes so they look about the same

pdf(snakemake@output[[7]], width = 11, height = 5)
# pdf("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/results/Figures/Fig2_MatingVsSNPs.pdf", width = 11, height = 5)
plot_grid(matingplot, NULL, chr5withline, ncol = 3, align = "h", rel_widths = c(2,0.1,2)) # This allows me to increase the space between plots
dev.off()

### Who are the hybrids in chromosome 5?
chr5pcaplot <- ggplot_build(chr5withline)
chr5pca <- chr5pcaplot$data[[1]]
hybridindex <- subset(chr5pca,x > -0.05 & x < 0) %>% rownames %>% as.numeric()
cat("The het-v haplotype recombinant samples are: ", samples[hybridindex])

# ============================
## Put them all together for Figure S4
# ============================

figureS4 <- plot_grid(PCAallgenome, PaPCA12, PaPCA13, PaPCA23, ncol = 1, labels = c('A', 'B', 'C', 'D'), rel_heights = c(0.6, 1,1,1), label_size = 20)

ggsave(snakemake@output[[8]], plot=figureS4, width = 11, height = 21)
# ggsave("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/8_SNPpop/results/Figures/FigS4_PCAs.pdf", plot=figureS4, width = 11, height = 21)

# ============================
print("DONE!", quote=FALSE)
