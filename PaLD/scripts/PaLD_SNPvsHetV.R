#!/usr/bin/env Rscript

### PaLD_SNPvsHetV: How linked are SNPs to the het-v locus along the genome?
#############################################################################
# Part of the Snakemake pipeline DiversityStats.smk
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2021-06-14
# =======================================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE)
library(vcfR)
library(cowplot)

# ============================
### Data
# ============================
# SNP data
vcf <- read.vcfR(snakemake@input$vcf, verbose = FALSE) # The vcf file

## Metadata
PopData <- read.csv(snakemake@input$PopData, header=TRUE, sep=",") # it only matters that it has het-v genotypes

# Output name
outputplot <- snakemake@output$LDplots

# IMPORTANT TO load aftr the output name, because it will replace it!!
heatmap <- load(snakemake@input$heatmap)
# The relevant plot is called `interchrs`. It didn't work in my local computer, but it did in Uppmax (where I created the object in the first place).


# ============================
### Functions
# ============================
# --- Calculate r2 between het genes ---
# Let's quantify the degree of gametic disequlibrium between the alelles of het-v and SNPs, using
# the square of the correlation of gene frequencies r2 (Hill & Robertson 1968).
# For this, I got inspired by the functions of the `genetics` package
print("Calculating LD between het-v and other het genes ...")
# Make a vector of presence absence of the allele of interest
haploallele.count <-function(gtps, allelo) {
  retval <- (gtps == allelo) %>% as.numeric() # Assume biallelic sites
  return(retval)
}

haploidLD <- function (g1, g2) {# Except g1 and g2 are not genetics::genotype() objects
  # Remove individuals with missing data
  bothgenes <- cbind.data.frame(g1,g2) %>% na.omit() %>% data.frame() # I had to use cbind.data.frame to avoid convertion of the factors to numeric
  # Calculate allele frequencies
  count.A <- bothgenes$g1 %>% table
  count.B <- bothgenes$g2 %>% table
  prop.A <- count.A/sum(count.A)
  prop.B <- count.B/sum(count.B)
  major.A <- names(prop.A)[which.max(prop.A)]
  major.B <- names(prop.B)[which.max(prop.B)]
  pA <- max(prop.A, na.rm = TRUE)
  pB <- max(prop.B, na.rm = TRUE)
  pa <- 1 - pA
  pb <- 1 - pB
  Dmin <- max(-pA * pB, -pa * pb) # for negative values of D (see Lewontin 1988)
  pmin <- pA * pB + Dmin
  Dmax <- min(pA * pb, pB * pa) # for positive values of D (see Lewontin 1988)
  pmax <- pA * pB + Dmax
  counts <- table(haploallele.count(g1, major.A), haploallele.count(g2, major.B))
  freqs <- counts/sum(counts)
  # The original LD.genotype function has a ML estimation equation because it 
  # doesn't know the phasing of the alleles in heterozygous. But this is haploid 
  # data, so we know. 
  # see https://cran.r-project.org/web/packages/genetics/genetics.pdf in page 25
  
  estD <- freqs[1,1]*freqs[2,2] - freqs[1,2]*freqs[2,1]
  if (estD > 0) estDp <- estD/Dmax else estDp <- estD/Dmin
  n <- sum(counts)
  corr <- estD/sqrt(pA * pB * pa * pb) # sensu Hill & Robertson (1968); this is r, so I need to get r2 by r^2
  # dchi <- (2 * n * estD^2)/(pA * pa * pB * pb) # see Lewontin 1988 for the Chi2 test of association
  # dpval <- 1 - pchisq(dchi, 1)
  
  # retval <- data.frame(D = estD, `D'` = estDp, 
  #                     r = corr, `R^2` = corr^2, n = n, `X^2` = dchi, `P-value` = dpval)
  # return(retval)
  return(corr^2)
  
  # # This totally depends on the genetics package
  # retval <- list(call = match.call(), D = estD, `D'` = estDp, 
  #                r = corr, `R^2` = corr^2, n = n, `X^2` = dchi, `P-value` = dpval)
  # class(retval) <- "LD"
  # retval
  #print(corr^2)
}

# ============================
### Make data frame of SNPs with their frequency
# ============================

# Remove biallelic SNPs
vcf <- vcf[is.biallelic(vcf),]

# Make a data frame
gt <- extract.gt(vcf, as.numeric = TRUE) # Extract genotypes is a matrix
nalt <- rowSums(gt) # How many individuals have the alternative allele?
nindiv <- dim(gt)[2] # How many individuals are in total? (notice there is no missing data to start)
vcffix <- getFIX(vcf) %>% data.frame # Extract the fixed part of the vcf

# Expand the data frame
allelecounts <- data.frame(CHROM = vcffix$CHROM, POS = vcffix$POS, snp = rownames(gt), nalt = nalt, maf = pmin(nalt, nindiv - nalt)/nindiv)
rownames(allelecounts) <- NULL

# Folded site frequency spectra
# ggplot(allelecounts %>% filter(maf > 0), aes(maf, fill = CHROM)) + geom_histogram() + facet_grid(CHROM ~ .)

# Let's remove low frequency alleles, say MAF > 0.02
cat("*** Filtering out SNPs with MAF < 0.02 ... ***\n")
cat("*** How many SNPs do we have?***\n")
allelecounts %>% dim # 145547      5
freqalleles <- allelecounts %>% filter(maf > 0.02)
cat("*** How many SNPs do we have after filtering?***\n")
dim(freqalleles) # 54486     5

## Make matrix of all SNPs for all strains
vcf_gt <- t(gt) %>% data.frame()

# Keep only sites with enough MAF
# common_vcf_gt <- vcf_gt %>% select(all_of(freqalleles$snp)) # in a higher version of dplyr I should have `select(all_of(freqalleles$snp))` because all_of is to prevent ambiguities with vectors (https://tidyselect.r-lib.org/reference/faq-external-vector.html)
common_vcf_gt <- vcf_gt[, names(vcf_gt) %in% freqalleles$snp] # More robust to use base r
common_vcf_gt <- cbind(Sample = rownames(common_vcf_gt), common_vcf_gt)
# rownames(vcf_gt) <- NULL

# Assign a het-v genotype to each individual
vcfall <- merge(common_vcf_gt, PopData %>% select(Sample, Het.v), by = "Sample")

### Calculate the LD between each SNP and het-v
cat("*** Calculating the LD between each SNP and het-v... ***\n")

# A little helper function
haploidLD_hetV<- function(col){
  retval <- haploidLD(col, vcfall$Het.v)
  return(retval)
}

# Calculate the R2 between each SNP and Het-v (it takes A WHILE)
R2perSNP <- apply(vcfall %>% select(-Sample), 2, haploidLD_hetV)

# Turn it into a data.frame
ldSNPs_all <- data.frame(R.2 = R2perSNP)
ldSNPs_all <- cbind(snp = rownames(ldSNPs_all), ldSNPs_all)
# Append Chromosome and positions
ldSNPs_all2 <- merge(freqalleles %>% select(CHROM, POS, snp), ldSNPs_all, by = "snp") %>% mutate(CHROM = factor(CHROM, levels = rev(c("chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7"))))

# How many SNPs per chromosome survived?
snpschrs <- ldSNPs_all2 %>% group_by(CHROM) %>% count

# Calculate media of each chromosome
# ldSNPs_all2 <- ldSNPs_all2 %>% group_by(CHROM) %>% mutate(medianR2 = median(R.2), q75 = quantile(R.2, prob=c(.75)))


LDsnpsHetV <- ggplot(ldSNPs_all2, aes(x = CHROM, y = R.2, fill = CHROM)) + geom_violin() +
  coord_flip() + # To avoid issues with "position_dodge requires non-overlapping x intervals" error, instead of just flipping x and y variables
  theme_bw() + theme(legend.position = "none", 
                     axis.title.y = element_blank(),
                     axis.text.y = element_text(size = 10, face = "bold")) +
  ylab(expression(LD~estimate~(r^2)~between~SNPs~and~the~italic("het-v")~locus)) +
  stat_summary(fun.y=median, geom='point', size=1, col='black') + # In higher dplyr versions, use `fun = 'median'` instead
  geom_text(data=snpschrs, aes(y = -0.06, x = CHROM, label=n), color="black", size = 3.2) + 
  scale_x_discrete(labels=c("chromosome_1" = "Chr 1", "chromosome_2" = "Chr 2", "chromosome_3" = "Chr 3",
                            "chromosome_4" = "Chr 4", "chromosome_5" = "Chr 5", "chromosome_6" = "Chr 6", "chromosome_7" = "Chr 7")) 


# Flipping x and y variables is easier and it works locally but not in Uppmax (probably due to older ggplot version)
# LDsnpsHetV <- ggplot(ldSNPs_all2, aes(x = R.2, y = CHROM, fill = CHROM)) + geom_violin() +
#   theme_bw() + theme(legend.position = "none", 
#                      axis.title.y = element_blank(),
#                      axis.text.y = element_text(size = 10, face = "bold")) +
#   xlab(expression(LD~estimate~(r^2)~between~SNPs~and~the~italic("het-v")~locus)) +
#   stat_summary(fun.y=median, geom='point', size=1, col='black') + # In higher dplyr versions, use `fun = 'median'` instead
#   geom_text(data=snpschrs, aes(x = -0.06, y = CHROM, label=n), color="black", size = 3.2) + 
#   scale_y_discrete(labels=c("chromosome_1" = "Chr 1", "chromosome_2" = "Chr 2", "chromosome_3" = "Chr 3",
#                             "chromosome_4" = "Chr 4", "chromosome_5" = "Chr 5", "chromosome_6" = "Chr 6", "chromosome_7" = "Chr 7"))
# 
cat("*** Putting plots together ***\n")
LDplots <- plot_grid(LDsnpsHetV, interchrs, labels = c('A', 'B'), ncol = 1, rel_heights = c(1, 2))

ggsave(outputplot, plot = LDplots, width = 19, height = 23, units = "cm")


cat("*** DONE! ***\n")
