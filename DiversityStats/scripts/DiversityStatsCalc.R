#!/usr/bin/env Rscript

### DiversityStatsCalc: Calculate diversity statistics of *Podospora anserina*
#############################################################################
# Part of the Snakemake pipeline DiversityStats.smk
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-07-12
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

# Repeats bed
badsites <- read.table(snakemake@input[[2]], header = TRUE, sep="\t")

outputname <- snakemake@output[[1]]
# -----
print("*** Reading data ***")
# Read entire vcf file and gff file
vcffolder <- dirname(snakemake@input[[3]]) # the location of the vcf
gfffolder <- dirname(snakemake@input[[4]]) # the location of the gff

Chrgc <- readData(vcffolder, format="VCF", gffpath = gfffolder, include.unknown = TRUE)

## Set different populations
Vclade <- PopData %>% filter(Het.v == "V") %>% .$Sample
Aclade <- PopData %>% filter(Het.v == "V1") %>% .$Sample

# poplist <- list(as.character(goodsamples), as.character(badsamples))
poplist <- list(as.character(Vclade), as.character(Aclade))

# ============================
# Define functions
# ============================

# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
#### Function to calculate a geometric mean, assuming that NAs are equivalent to 0's 
# (in my meanR2 function it doesn't matter because the NAs will be catched first)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#### Function to calculate the mean r^2 value of a window in a GENOME.class object (transformed into windows)
meanR2 <- function(GENOME.class.slide, popindex = 1){
  R2wind <- c()
  for (region in 1:length(GENOME.class.slide@region.names)){
    # print(region)
    # print(GENOME.class.slide@SLIDE.POS[[region]])

    # Sometimes the region doesn't have any SNPs and no LD is calculated; in those cases, there is not 
    # even a second population so popindex can be out of bounds.
    regionmatrix <- tryCatch(
      {
        # [[x]][[y]] x:region, y:population 
        regionmatrix <- GENOME.class.slide@region.stats@linkage.disequilibrium[[region]][[popindex]]
        
        if(length(regionmatrix) == 1){ # Sometimes the result is just NaN
          if(is.na(regionmatrix)) {regionmatrix <- NULL}
        }
        
        regionmatrix # this is the output
        # return(regionmatrix)
      },
      error = function(cond){  # In case there is no second (or third ...) pop
        # message(cond) # Print the error message
        regionmatrix <- NULL
        return(regionmatrix)
      }
    )
    
    ## Calculate the mean of that region
    if(is.null(regionmatrix)){ # There were no SNPs in the window, probably
      R2wind <- c(R2wind, NA)
    } else{
      regionr2 <- regionmatrix %>% t() %>% data.frame
      ## Arithmetic mean
      # meanr2 <- suppressWarnings(mean(regionr2$R2)) # The NAs will give a warning everytime
      ## Geometric mean
      meanr2 <- gm_mean(regionr2$R2)
      R2wind <- c(R2wind, meanr2)
    }
  }
  R2wind <- as.numeric(R2wind)
  return(R2wind)
}

#### Function to calculate the amount of missing data per window
## The idea is to reduce the size of each window with the overlapping regions in badsites

winmissing <- function(GENOME.class.slide, badsites){ #, minlen = 5000
  # Get coordinates of the windows for this GENOME.class.slide object
  coords <- c()
  for (win in GENOME.class.slide@region.names) {
    start <- strsplit(win," ")[[1]][c(1)] %>% as.numeric()
    end <- strsplit(win," ")[[1]][c(3)] %>% as.numeric()
    coords <- rbind(coords, cbind(start, end))
  }
  coords <- data.frame(coords)
  
  winlens <- c()
  for (win in seq(1, nrow(coords))){
    start1 <- coords[win,1] %>% as.numeric()
    end1 <- coords[win,2] %>% as.numeric()
    
    fulllenwin <- end1 - start1 + 1 # Plus one because the range includes the ends
    for (bad in seq(1, nrow(badsites))) {
      start2 <- badsites[bad,2] %>% as.numeric()
      end2 <- badsites[bad,3] %>% as.numeric()
      
      # Is there overlap between ranges?
      if ((start1 <= end2) & (end1 >= start2)) { # Yes
        # How much?
        if (start2 <= start1 & end2 <= end1) { # The bad range starts before the window and ends within
          fulllenwin <- fulllenwin - (end2 - start1 + 1)
        } else if (start2 < start1 & end2 > end1) { # The bad range is larger and includes the window
          fulllenwin <- 0
          break
        } else if (start2 > start1 & end2 < end1) { # The bad range is fully contained in the window
          fulllenwin <- fulllenwin - (end2 - start2 + 1)
        } else if (start2 > start1 & end2 > end1) { # The bad range starts within the window but ends after
          fulllenwin <- fulllenwin - (end1 - start2 + 1)
        }
      } else if (start2 > end1){ # Are the bad ranges after the window? then there is no point in continue
        break
      }
    }
    
    ## I commented it because I deal with the allowed number of sizes later
    # # Only windows of a minimum size are accepted
    # if (fulllenwin <= minlen){
    #   fulllenwin <- 0
    # }
    
    winlens <- c(winlens, fulllenwin)
  }
  return(winlens)
}

### Function to create a data frame of pop gen statistics for an input GENOME.class object with one population
statschrs <- function(GENOME.class, badsites, width= 10000, jump = 1000, pop1 = "All", r2 = FALSE, indvstats = TRUE, precalcwinlens = NULL){
  # pop1 is just a name for the population
  # badsites is a dataframe (bed-like) with ranges in the reference genome that should be set as missing data
  
  cat('\n')
  cat("*** Window transformation ***\n")
  GENOME.class.slide <- sliding.window.transform(GENOME.class, width = width, jump = jump, type=2) # 1 scan only biallelic positions (SNPs), 2 scan the genome
  
  # ------- Individual stats
  
  # The slot GENOME.class@region.names will store the genomic regions of each window
  # as a character string. To convert those strings into a genomic numeric
  # position we can apply the following
  genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    val <- mean(as.numeric(split))
    return(val)
  })
  
  # Since it takes so long to calculate, it's nice to have the option to input the winlens data frame from the beginning 
  if(is.null(precalcwinlens)){
    cat("*** Calculating missing data per window ***\n")
    winlens <- winmissing(GENOME.class.slide, badsites) #, minlen = width/2)
  } else{
    winlens <- precalcwinlens
  }
  
  if(indvstats!=FALSE){
    # Tajima
    cat("*** Tajima's D ***\n")
    GENOME.class.slide <- neutrality.stats(GENOME.class.slide, FAST = TRUE) # FAST is set to TRUE. This will speed up calculations but might be a bit unstable in some cases.
    
    # Diversity stats
    cat("*** Nucleotide diversity ***\n")
    GENOME.class.slide <- diversity.stats(GENOME.class.slide)
    nucdiv <- GENOME.class.slide@nuc.diversity.within # Same as GENOME.class.slide@theta_Tajima
    nucdiv <- nucdiv/winlens # the values have to be normalized by the number of nucleotides in each window
    
    # Get the raw Pi
    Pi <- data.frame(position = genome.pos, pop1 = nucdiv[,1]) %>% melt(id = c("position"))
    names(Pi) <- c("position", "Pop", "Pi")
    Pi$Pi[which(Pi$Pi == Inf)] <- NA # Remove infinite values
    
    # Get Waterson's theta
    watersontheta <- GENOME.class.slide@theta_Watterson/winlens
    waterson <- data.frame(position = genome.pos, pop1 = watersontheta[,1]) %>% melt(id = c("position"))
    waterson$value[which( waterson$value == Inf)] <- NA # Remove infinite values
    
    # Get the Tajima's D
    TD <- data.frame(position = genome.pos, pop1 = GENOME.class.slide@Tajima.D[,1]) %>% melt(id = c("position"))
    
    # Put them together
    windsdf <- cbind(Pi, Tajima = TD$value, Theta_W = waterson$value) 
    # Melt them to work on ggplot better
    windsdf <- windsdf %>% melt(id = c("Pop", "position"))
    
    # Rename the populations
    windsdf$Pop <- pop1
    levels(windsdf$Pop) <- c(pop1)
  } else {
    windsdf <- data.frame(Pop = rep(pop1, length(genome.pos)), position = genome.pos, row.names = NULL)
  }
  
  # Put them together
  # pairwisestats <- data.frame(Dxy = dxypop$pop1.pop2, Fst = pairwise.FST[,1])
  
  # ------- 
  ## Linkage disequilibrium
  r2df <- c(NULL, NULL)
  if(r2!=FALSE){
    cat("*** Linkage disequilibrium ***\n")
    # To get r2 (Hill & Robertson, 1968)
    GENOME.class.slide <- calc.R2(GENOME.class.slide, lower.bound = 0.05) # TAKES FOREVER, like an hour for chr4 without a lower.bound!!!
    # GENOME.class.slide <- calc.R2(GENOME.class.slide, lower.bound = 3/(get.individuals(GENOME.class)[[1]] %>% length())) # TAKES FOREVER, like an hour for chr4 without a lower.bound!!!
    # Result 
    # [[x]][[y]] x:region, y:population 
    # GENOME.class.slide@region.stats@linkage.disequilibrium[[1]][[1]] %>% t() %>% head
    
    # Get the mean r2 for each window
    pop1r2 <- meanR2(GENOME.class.slide, popindex = 1) 
    # Make a data grame
    pop1r2df <- data.frame(Pop = rep(pop1, length(genome.pos)), position = genome.pos, variable = "r2", value = pop1r2, row.names = NULL)
    # attach to the main data frame
    windsdf <- rbind(windsdf, pop1r2df)
  }
  
  # Add a column with the number of effective sites pair window
  windsdf <- cbind(windsdf, winlens)
  
  # ------- Reporting
  # Results:
  # * The modified GENOME.class with windows transformation
  # * The dataframe with the results of all samples together
  results <- list(GENOME.class.slide, windsdf)
  return(results)
}

### Function to create a data frame of pop gen statistics for an input GENOME.class object using windows with two populations
statschrs2pops <- function(GENOME.class, poplist, badsites, width= 10000, jump = 1000, pop1 = "V", pop2 = "V1", r2 = FALSE, indvstats = TRUE, precalcwinlens = NULL){
  # pop1 and pop2 are just names for the populations
  # badsites is a dataframe (bed-like) with ranges in the reference genome that should be set as missing data
  cat("*** Setting populations ***\n")
  GENOME.class <- set.populations(GENOME.class, poplist, diploid = FALSE)
  
  cat('\n')
  cat("*** Window transformation ***\n")
  GENOME.class.slide <- sliding.window.transform(GENOME.class, width = width, jump = jump, type=2) # 1 scan only biallelic positions (SNPs), 2 scan the genome
  
  # ------- Individual stats
  
  # The slot GENOME.class@region.names will store the genomic regions of each window
  # as a character string. To convert those strings into a genomic numeric
  # position we can apply the following
  genome.pos <- sapply(GENOME.class.slide@region.names, function(x){
    split <- strsplit(x," ")[[1]][c(1,3)]
    val <- mean(as.numeric(split))
    return(val)
  })
  
  # Since it takes so long to calculate, it's nice to have the option to input the winlens data frame from the beginning 
  if(is.null(precalcwinlens)){
    cat("*** Calculating missing data per window ***\n")
    winlens <- winmissing(GENOME.class.slide, badsites) #, minlen = width/2)
  } else{
    winlens <- precalcwinlens
  }
  
  if(indvstats!=FALSE){
    # Tajima
    cat("*** Tajima's D ***\n")
    GENOME.class.slide <- neutrality.stats(GENOME.class.slide, FAST = TRUE) # FAST is set to TRUE. This will speed up calculations but might be a bit unstable in some cases.
    
    # Diversity stats
    cat("*** Nucleotide diversity ***\n")
    GENOME.class.slide <- diversity.stats(GENOME.class.slide)
    nucdiv <- GENOME.class.slide@nuc.diversity.within # Same as GENOME.class.slide@theta_Tajima
    nucdiv <- nucdiv/winlens # the values have to be normalized by the number of nucleotides in each window
    
    # Get the raw Pi
    Pi <- data.frame(position = genome.pos, pop1 = nucdiv[,1], pop2 = nucdiv[,2]) %>% melt(id = c("position"))
    names(Pi) <- c("position", "Pop", "Pi")
    Pi$Pi[which(Pi$Pi == Inf)] <- NA # Remove infinite values
    
    # Get Waterson's theta
    watersontheta <- GENOME.class.slide@theta_Watterson/winlens
    waterson <- data.frame(position = genome.pos, pop1 = watersontheta[,1], pop2 = watersontheta[,2]) %>% melt(id = c("position"))
    waterson$value[which( waterson$value == Inf)] <- NA # Remove infinite values
    
    # Get the Tajima's D
    TD <- data.frame(position = genome.pos, pop1 = GENOME.class.slide@Tajima.D[,1], pop2 = GENOME.class.slide@Tajima.D[,2]) %>% melt(id = c("position"))
    
    # Put them together
    windsdf <- cbind(Pi, Tajima = TD$value, Theta_W = waterson$value) 
    # Melt them to work on ggplot better
    windsdf <- windsdf %>% melt(id = c("Pop", "position"))
    
    # Rename the populations
    windsdf$Pop <- windsdf$Pop %>% gsub("pop1", pop1, .)
    windsdf$Pop <- windsdf$Pop %>% gsub("pop2", pop2, .)
    levels(windsdf$Pop) <- c(pop1, pop2)
  } else {
    windsdf <- data.frame(Pop = rep(pop1, length(genome.pos)), position = genome.pos, row.names = NULL)
  }
  
  # ------- Stats between populations
  cat("*** Dxy ***\n")
  GENOME.class.slide <- diversity.stats.between(GENOME.class.slide)
  dxypop <- data.frame(GENOME.class.slide@nuc.diversity.between/winlens) # Between pops
  dxypop$pop1.pop2[which(dxypop$pop1.pop2 == Inf)] <- NA
  # dxypop <- data.frame(t(GENOME.class.slide@nuc.diversity.between)/width) # Between pops
  
  cat("*** Fst ***\n")
  GENOME.class.slide <- F_ST.stats(GENOME.class.slide, mode="nucleotide")
  # Get the pairwise nucleotide FST for the two populations
  pairwise.FST <- data.frame(t(GENOME.class.slide@nuc.F_ST.pairwise))
  # Remove negative values, maybe they are there because estimators of Hs and Ht are used? (like in Gst_Hedrick {mmod})
  pairwise.FST$pop1.pop2[which(pairwise.FST$pop1.pop2 < 0)] <- 0
  
  # Put them together
  # pairwisestats <- data.frame(Dxy = dxypop$pop1.pop2, Fst = pairwise.FST[,1])
  pairwisestats <- data.frame(Pop = c(rep("VvsV1", length(genome.pos))), position = genome.pos, variable = c(rep("Fst", length(genome.pos)), rep("Dxy", length(genome.pos))), value = c(pairwise.FST[,1], dxypop$pop1.pop2), row.names = NULL)
  
  # ------- 
  ## Linkage disequilibrium
  r2df <- c(NULL, NULL)
  if(r2!=FALSE){
    cat("*** Linkage disequilibrium ***\n")
    # To get r2 (Hill & Robertson, 1968)
    GENOME.class.slide <- calc.R2(GENOME.class.slide, lower.bound = 0.05) # TAKES FOREVER, like an hour for chr4 without a lower.bound!!!
    # GENOME.class.slide <- calc.R2(GENOME.class.slide, lower.bound = 3/(get.individuals(GENOME.class)[[1]] %>% length())) # TAKES FOREVER, like an hour for chr4 without a lower.bound!!!
    # Result 
    # [[x]][[y]] x:region, y:population 
    # GENOME.class.slide@region.stats@linkage.disequilibrium[[1]][[1]] %>% t() %>% head
    
    # Get the mean r2 for each window
    pop1r2 <- meanR2(GENOME.class.slide, popindex = 1)
    pop2r2 <- meanR2(GENOME.class.slide, popindex = 2)
    
    # Make a data grame
    r2df <- data.frame(Pop = c(rep(pop1, length(genome.pos)), rep(pop2, length(genome.pos))), position = genome.pos, variable = "r2", value = c(pop1r2, pop2r2), row.names = NULL)
    # attach to the main data frame
    windsdf <- rbind(windsdf, r2df)
    
  }
  
  # ------- 
  cat("### Putting stats together ... ###\n")
  windsdf <- rbind(windsdf, pairwisestats) 
  
  # Add a column with the number of effective sites pair window
  windsdf <- cbind(windsdf, winlens)
  
  # ------- Reporting
  # Results:
  # * The modified GENOME.class with windows transformation
  # * The dataframe with the results within the two populations samples together
  # * The dataframe with the results of comparing two pops
  # results <- list(GENOME.class.slide, windsdf, pairwisestats)
  results <- list(GENOME.class.slide, windsdf)
  return(results)
}

### Function to plot several statistic along the chromosome, assuming only two populations, in a pre-processed dataframe from statschrs() or statschrs2pops()
# Notice it will return a plot in facet_grid format
statschrs.plotter <- function(windsdf, title = "Chromosome 4"){ 
  p <- ggplot(windsdf, aes(x = position, y = value, fill = Pop, colour = Pop)) + 
    geom_line(alpha = 0.9) + 
    facet_grid(variable ~., scales="free_y") +
    labs(title = title, x = "Position (bp)") + 
    theme_bw() + # Remove background color
    scale_x_continuous(expand = c(0.02, 0.02)) + # Remove the white space on the sides
    theme(plot.title = element_text(size = rel(1.2), face="bold", family = "Helvetica", hjust = 0.5), axis.title.y=element_blank()) # face="bold.italic" # , panel.grid = element_blank(), panel.border = element_blank())
  return(p)
}

### Function to calculate all the diversity stats, within and between populations and from all samples
# It returns a single data frame
statswithinandbetween <- function(GENOME.class, poplist, badsites, width = 10000, jump = 1000, r2 = FALSE){
  cat("*** Slicing genome class object ***\n")
  GENOME.class.slide <- sliding.window.transform(GENOME.class, width = width, jump = jump, type=2) # 1 scan only biallelic positions (SNPs), 2 scan the genome
  cat("*** Calculating missing data per window ***\n")
  winlens <- winmissing(GENOME.class.slide, badsites)
  
  cat("*** Calculating stats for all samples ***\n")
  Chrgc.sc.all <-  statschrs(GENOME.class, badsites, width = width, jump = jump, r2 = r2, precalcwinlens = winlens) 
  Chrgc.sc.all.df <- Chrgc.sc.all[[2]]
  
  cat("*** Calculating stats within and between pops ***\n")
  Chrgc.sc <-  statschrs2pops(GENOME.class, poplist, badsites, width = width, jump = jump, pop1 = "V", pop2 = "V1", r2 = r2, precalcwinlens = winlens) 
  Chrgc.sc.df <- Chrgc.sc[[2]]
  
  result <- rbind(Chrgc.sc.all.df, Chrgc.sc.df)
  return(result)
}

# ============================
# Analysis
# ============================
allchrstats <- statswithinandbetween(Chrgc, poplist, badsites, width = 10000, jump = 1000, r2 = FALSE) # I decided to turn off the linkage part
write.csv(allchrstats, file = outputname, row.names = FALSE)

# # ============================
# ### For some reason this crashes for chromosomes 2 and 7 in UPPMAX. If I run it locally in my computer it works fine...
# # ============================