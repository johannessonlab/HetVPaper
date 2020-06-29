#!/usr/bin/env Rscript

# DiversityStats_vcfR_plotter: Plot the coverage distribution of the *P. anserina* samples
#############################################################################
# Part of the Snakemake pipeline DiversityStats.smk 
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-07-17
# Version 1
#############################################################################

# ============================
# Check input
# ============================
vcffile <- snakemake@input[[1]]
plotcov <- snakemake@output[[1]]
tablecov <- snakemake@output[[2]]

lquantile = snakemake@params$lquantile
uquantile = snakemake@params$uquantile

# ============================
# Load the necessary libraries
# ============================
library(vcfR)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE) # So it won't print unnecessary stuff

# ============================
# Reading the data
# ============================
vcf <- read.vcfR(vcffile, verbose = FALSE)

# ============================
# Prepare the functions
# ============================
# vignette('sequence_coverage')
coolviolins <- function(dp){
  if( require(reshape2) & require(ggplot2) ){
    dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
    dpf <- dpf[ dpf$Depth > 0,]
    p <- ggplot(dpf, aes(x=Sample, y=Depth)) + geom_violin(fill="#C0C0C0", adjust=1.0,
                                                           scale = "count", trim=TRUE)
    p <- p + theme_bw()
    p <- p + theme(axis.title.x = element_blank(), 
                   axis.text.x = element_text(angle = 60, hjust = 1))
    p <- p + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black")
    p <- p + scale_y_continuous(trans=scales::log2_trans(), breaks=c(1, 10, 100, 800))
    p
    return(p)
  } else {
    message("The packages reshape2 and ggplot2 are required for this example but do not appear
          to be installed. Please use install.packages(c('reshape2', 'ggplot2', 'scales')) if you would
          like to install them.")
  }
}

# ============================
# Analysis
# ============================
## Print violin plots of coverage
dp <- extract.gt(vcf, element="DP", as.numeric = TRUE) # Make data frame of depth of coverage

ggsave(plotcov, plot = coolviolins(dp), width = 30, height = 10, units = "cm")

# Print a summary table of the coverage distribution per sample
# I chose 0.985 out of trial and error. I tested 0.9-0.99
dpsum <- data.frame(sample = colnames(dp), median = apply(dp, 2, median, na.rm = TRUE), mean = apply(dp, 2, mean, na.rm = TRUE), lowquantile = apply(dp, 2, quantile, probs = lquantile, na.rm = TRUE), upquantile = apply(dp, 2, quantile, probs = uquantile, na.rm = TRUE))

write.table(dpsum, file = tablecov, sep = "\t", row.names=FALSE, quote=FALSE)

