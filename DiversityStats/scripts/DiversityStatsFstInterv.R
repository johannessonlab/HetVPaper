#!/usr/bin/env Rscript

### DiversityStatsFstInterv: Calculating confidence intervals for Fst in *Podospora anserina* chromosomes
#############################################################################
# Part of the Snakemake pipeline DiversityStats.smk
# Based on PaFstConfintervals.R
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-08-01
# Version 1
# =======================================
# https://github.com/tonig-evo/workshop-popgenome
# https://www.rdocumentation.org/packages/PopGenome/versions/2.2.4/topics/neutrality.stats-methods
# http://popgenome.weebly.com/ 
# https://www.biostars.org/p/121428/

library(PopGenome)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(reshape2)

# ============================
# Reading the data
# ============================
# Metadata
PopData <- read.csv(snakemake@input[[1]], header = TRUE, sep=",")
# PopData <- read.csv("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/9b_DiversityStatsWa/data/SNPpopMetadata.csv", header = TRUE, sep=",")
# PopData <- read.csv("/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/9b_DiversityStatsWa/data/PaAllelesPaper.csv", header = TRUE, na.strings = "NA", sep=",")

# -----
print("*** Reading data ***")
# Read entire vcf file and gff file
vcffolder <- dirname(snakemake@input[[2]]) # the location of the vcf
# vcffolder <- "/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/9_DiversityStats/data/chr2" # the location of the vcf
gfffolder <- dirname(snakemake@input[[3]]) # the location of the gff
# gfffolder <- "/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/9_DiversityStats/data/gffs" # the location of the gff

Chrgc <- readData(vcffolder, format="VCF", gffpath = gfffolder, include.unknown = TRUE)

## Set different populations
Vclade <- PopData %>% filter(Het.v == "V") %>% .$Sample
Aclade <- PopData %>% filter(Het.v == "V1") %>% .$Sample

# poplist <- list(as.character(goodsamples), as.character(badsamples))
poplist <- list(as.character(Vclade), as.character(Aclade))

# CONSTANTS
Nruns <- 1000

# ============================
# Define functions
# ============================

#### Function to create random subsets of samples WITHOUT? replacement
randomgroups <- function(samplesvect, n1 = poplist[[1]] %>% length()){ # Use the size of the real mating groups
  # Make fake subpops
  sample1 <- as.character(samplesvect) %>% sample(size = n1)
  sample2 <- as.character(samplesvect[!(as.character(samplesvect) %in% sample1)])
  # Make them a list
  newpoplist <- list(sample1, sample2)
  return(newpoplist)
}

#### Function to produce Fst values for an input population list; the nsim argument is there to control the name of the output column
FstChr <- function(GENOME.class, poplist = randomgroups(PopData$Sample), width = 10000, jump = 1000, nsim = 1){
  cat("*** Setting populations ***\n")
  # Re-difine pops based on the random assignment of samples to two groups
  GENOME.class <- set.populations(GENOME.class, poplist, diploid = FALSE)
  
  cat('\n')
  cat("*** Window transformation ***\n")
  # Run again the window transformation to match the pops
  GENOME.class.slide <- sliding.window.transform(GENOME.class, width = width, jump = jump, type=2) # 1 scan only biallelic positions (SNPs), 2 scan the genome
  
  cat("*** Fst ***\n")
  GENOME.class.slide <- F_ST.stats(GENOME.class.slide, mode="nucleotide")
  # Get the pairwise nucleotide FST for the two populations
  pairwise.FST <- data.frame(t(GENOME.class.slide@nuc.F_ST.pairwise))
  
  # Remove negative values, maybe they are there because estimators of Hs and Ht are used? (like in Gst_Hedrick {mmod})
  pairwise.FST$pop1.pop2[which(pairwise.FST$pop1.pop2 < 0)] <- 0
  names(pairwise.FST) <- paste0("Replica", nsim, sep = "")
  return(pairwise.FST)
}
# FstChr(GENOME.class) %>% head

### Function to calculate confidence intervals of Fst based on a permutation test; dependant on FstChr(), randomgroups() and goodsamples/clades
FstCIs <- function(GENOME.class, repfile, Nruns = 100, width = 10000, jump = 1000){
  ## Define the regions once without the populations set yet
  GENOME.class.slide <- sliding.window.transform(GENOME.class, width = width, jump = jump, type=2) # 1 scan only biallelic positions (SNPs), 2 scan the genome
  
  ## Get the mid positions of each interval
  # The slot GENOME.class@region.names will store the genomic regions of each window
  # as a character string. To convert those strings into a genomic numeric
  # position we can apply the following
  
  genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    val <- mean(as.numeric(split))
    return(val)
  })
  
  starts <- sapply(GENOME.class.slide@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    start <- split[1]
    return(start)
  })
  
  ends <- sapply(GENOME.class.slide@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    end <- split[2]
    return(end)
  })
  
  # Make the permutation test
  # [Rows sites, Columns Simulation]
  FstRounds <- matrix(nrow = length(genome.pos), ncol = Nruns, dimnames = list(genome.pos, 1:Nruns))
  for (i in 1:Nruns){
    cat("*** Round", i, "***\n")
    FstRun <- FstChr(GENOME.class, poplist = randomgroups(PopData$Sample), width = 10000, jump = 1000, nsim = i)
    FstRounds[,i] <- FstRun %>% as.matrix
  }
  
  # Save the replicas elsewhere
  write.csv(cbind(starts, ends, genome.pos, FstRounds), file = repfile, row.names = FALSE)
  
  cat("*** Calculating observed values ***\n")
  # Make a data frame with the observed values
  ActualFst <- data.frame(start = as.vector(starts), end = as.vector(ends), position = genome.pos, variable = rep("Fst", length(genome.pos)),
                          value = FstChr(GENOME.class, poplist = poplist, width = width, jump = jump), row.names = NULL)
  names(ActualFst)[5] <- "value" # Re-name because it didn't work above
  
  
  # Calculate confidence intervals for each site (window)
  Q1Fst <- c()
  Q3Fst <- c()
  maxFst <- c()
  for (i in 1:length(genome.pos)){
    # Calculate confidence intervals
    quantiles <- quantile(FstRounds[i,], probs = c(0.025, 0.975), na.rm = TRUE, type = 7) # Two tails (default quantile type)
    Q1Fst <- c(Q1Fst, quantiles[1])
    Q3Fst <- c(Q3Fst, quantiles[2])
    maxFst <- c(maxFst, max(FstRounds[i,]))
  }
  
  # Put all together
  FstData <- cbind(ActualFst, Q1Fst, Q3Fst, maxFst)
  
  return(FstData)
}

# ============================
# Analysis
# ============================
# Chrgc.Fst <- FstCIs(Chrgc, "/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/9_DiversityStats/data/FstReplicas/chr2.csv", Nruns = Nruns)
Chrgc.Fst <- FstCIs(Chrgc, snakemake@output[[1]], Nruns = Nruns)
# write.csv(Chrgc.Fst, file = "FstCIs/FstRandom_Chr4.csv", row.names = FALSE)
write.csv(Chrgc.Fst, file = snakemake@output[[2]], row.names = FALSE) # write the digested data.frame
