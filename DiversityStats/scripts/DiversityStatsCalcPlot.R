#!/usr/bin/env Rscript

### DiversityStatsCalcPlot: Plot the diversity statistics of *Podospora anserina* calculated in DiversityStatsCalc.R
#############################################################################
# Part of the Snakemake pipeline DiversityStats.smk
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-07-15
# Version 2
# =======================================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE)
library(reshape2)
library(tidyverse) # To rearrange factors easily
library(cowplot) # For plot_grid
library(stringr) # To print legend title with new lines
library(grid) # for textGrob()
library(gridExtra) # for grid.arrange()

# ============================
# Reading the data
# ============================
### SNP data in vcf files
PopDataChr1 <- read.csv(snakemake@input[[1]], header = TRUE, sep=",")
PopDataChr2 <- read.csv(snakemake@input[[2]], header = TRUE, sep=",")
PopDataChr3 <- read.csv(snakemake@input[[3]], header = TRUE, sep=",")
PopDataChr4 <- read.csv(snakemake@input[[4]], header = TRUE, sep=",")
PopDataChr5 <- read.csv(snakemake@input[[5]], header = TRUE, sep=",")
PopDataChr6 <- read.csv(snakemake@input[[6]], header = TRUE, sep=",")
PopDataChr7 <- read.csv(snakemake@input[[7]], header = TRUE, sep=",")

### Fst calculations with the confidence intervals
FstConfIntvChr1 <- read.csv(snakemake@input[[8]], header = TRUE, sep=",")
FstConfIntvChr2 <- read.csv(snakemake@input[[9]], header = TRUE, sep=",")
FstConfIntvChr3 <- read.csv(snakemake@input[[10]], header = TRUE, sep=",")
FstConfIntvChr4 <- read.csv(snakemake@input[[11]], header = TRUE, sep=",")
FstConfIntvChr5 <- read.csv(snakemake@input[[12]], header = TRUE, sep=",")
FstConfIntvChr6 <- read.csv(snakemake@input[[13]], header = TRUE, sep=",")
FstConfIntvChr7 <- read.csv(snakemake@input[[14]], header = TRUE, sep=",")

## The het genes historical data
hetrecord <- read.csv(snakemake@input[[15]], header = TRUE, na.strings = "NA", sep=",")
hetrecord$Het.Z <- as.factor(hetrecord$Het.Z)

# Rearrange for the genes
alleledata <- melt(hetrecord %>% select(SampleID, Year, Het.v, Het.Z, Het.r_phen, Het.c, Het.c_phen, Het.s, SimpleID), id = c("SampleID", "Year"))
names(alleledata)[3:4] <- c("Het", "Allele")

PopDataChr5noHybrids <- read.csv(snakemake@input[[16]], header = TRUE, sep=",")
# gff_file <- snakemake@input[[17]]
# TE_file <- snakemake@input[[18]]

# ------

## The coordinates of the genes, used by the plotting functions
genesall <- data.frame(locus = c("Spok2", "MAT", rep("Centromere", 7), 
                                 "het-z", "het-r", "het-d", "het-s", "het-c", "nwd1", # Deskalov et al 2012: het-s is next to NWD2
                                 "hnwd1", "het-e", "mod-A", "het-V", "nwd3", "idi-2",
                                 "hnwd3", "Spok3/4", "Spok3/4", "Spok3/4"),
                       position = c(1094363, 7345878.5, 4479205, 236468.5, 675115.5, 1236388.5, 2062808.5, 238150, 3562141.5,
                                    3957511.5, 2683360.5, 3056777, 208504, 495416,  3983565.5,
                                    436378.5, 691354.5, 1568320, 1294132.5, 1392335, 1641319,
                                    712318.5, 355881, 3328395.5, 900387.5),
                       chr = c("chromosome_5", "chromosome_1","chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7",
                               "chromosome_1", "chromosome_2", "chromosome_2", "chromosome_3", "chromosome_3", "chromosome_3",
                               "chromosome_4", "chromosome_4", "chromosome_4", "chromosome_5", "chromosome_5", "chromosome_5",
                               "chromosome_7", "chromosome_3", "chromosome_5", "chromosome_5"),
                       infig = c(1,1,1,1,1,1,1,1,1,  # Use this to filter out loci
                                 1,1,1,1,1,0,
                                 1,1,0,1,0,0,
                                 1,1,1,1),
                       shapes = c(10, 15, rep(19,7), 
                                  17, 24, 4, 5, 6, 8,
                                  12, 18, 13, 25, 2, 11,
                                  9, 1, 1, 1),
                       value = 0)

genes <- genesall %>% filter(infig == 1)

# Manual
maxPi <- 0.02

# ============================
# Functions
# ============================
# To print the scale with more digits
scaleFUN <- function(x) sprintf("%.3f", x)

# Plot Pi, Waterson's Theta, and Tajima's D, but allowing for geom_area() in one and not the other
allpitajima <- function(PopData, thischr = "chromosome_5", thistitle = "Chromosome 5", minlen = 5000, legend = "legend"){
  PopDataFiltered <- PopData %>% filter(Pop == "All") %>% filter(variable != "r2") 
  PopDataFiltered$value[which(PopDataFiltered$winlens <= minlen )] <- NA 
  
  piplot <- ggplot(PopDataFiltered %>% filter(variable != "Tajima"), aes(x = position, y = value, colour = variable)) + 
    theme_bw() + # Remove background color
    labs(title = thistitle, y = expression(pi*" or "*theta["W"])) +
    scale_x_continuous(expand = c(0.02, 0.02)) + # Remove the white space on the sides
    theme(plot.title = element_text(size = rel(1.5), face="bold", family = "Helvetica", hjust = 0.5), 
          # plot.margin = margin(t = 0, r = 10, b = 0, l = 3, unit = "pt"),
          axis.title.y= element_text(size = rel(1.5),),
          axis.title.x=element_blank(), 
          axis.text = element_text(size = rel(1)),
          legend.key=element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(size=rel(0.9)),
          legend.position=c(0.94, 0.56),
          legend.spacing.y = unit(0.005, 'cm'), # Put the items of the legend closer to each other
          legend.background=element_blank()) +
    scale_color_manual(values= c("firebrick1", "blue4"), labels=c(expression(pi["  "]), expression(theta["W"])) , guide=legend) + # The space in pi is to center it better...
    ylim(0, round(maxPi,digits=2)) +
    geom_line(alpha = 0.7) + 
    geom_point(data = genes %>% filter(chr == thischr) %>% filter(locus != "nwd3"), aes(x = position, shape = locus, fill = locus), colour = "black", size = 3) + guides(fill=FALSE) + # %>% filter(locus != "het-V") 
    # scale_shape_manual(values=c(19, 2, 5, 1, 17, 23))
    scale_shape_manual(values= genes[order(genes$locus),] %>% filter(chr == thischr) %>% filter(locus != "nwd3") %>% .$shapes) # Give it specific shapes (but re-order the dataframe so it maches) #  %>% filter(locus != "het-V")
  
  tajimaplot <- ggplot(PopDataFiltered %>% filter(variable == "Tajima"), aes(x = position, y = value, colour = variable, fill = variable)) + 
    geom_line(alpha = 0.4) + 
    theme_bw() + # Remove background color
    labs(x = "Position (bp)", y = "Tajima's D") + 
    scale_x_continuous(expand = c(0.02, 0.02)) + # Remove the white space on the sides
    theme(plot.title = element_text(size = rel(1.5), face="bold", family = "Helvetica", hjust = 0.5), 
          # plot.margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"),
          axis.title.y= element_text(size = rel(1.5),),
          axis.title.x= element_text(size = rel(1.3)),
          axis.text = element_text(size = rel(1)),
          legend.key=element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(size=rel(0.9)),
          legend.position=c(0.94, 0.9),  # xy
          legend.background=element_blank()) +
    scale_color_manual(values="darkorange3", labels= "Tajima's D", guide=legend) +
    geom_area(alpha = 0.5, fill = "darkorange3") +
    scale_y_continuous(labels=scaleFUN, limits = c(-3, 5)) +
    geom_point(data = genes %>% filter(chr == thischr) %>% filter(locus != "nwd3"), aes(x = position, shape = locus, fill = locus), colour = "black", size = 3,  # %>% filter(locus != "het-V")
               show.legend = FALSE) + guides(fill=FALSE) + # The guides remove part of the legend
    scale_shape_manual(values= genes[order(genes$locus),] %>% filter(chr == thischr) %>% filter(locus != "nwd3") %>% .$shapes) # Give it specific shapes (but re-order the dataframe so it maches) # %>% filter(locus != "het-V")
  
  together <- plot_grid(piplot, tajimaplot, ncol = 1, align = "v") # the align part ensures the x axis matches
  return(together)
}

# A function to calculate the allele frequencies of a given het gene data frame
yearlycounts <- function(hetdf, minsample = 5){
  # Make empty dataframe
  PopAlleleCount <- data.frame(year = NA, Allele = NA, Freq = NA, Proportion = NA)
  
  # How many alleles are there in total?
  allAlleles <- hetdf$Allele %>% na.omit() %>% unique() 
  
  # Loop trough the years
  for (year in hetdf$Year %>% na.omit %>% unique() %>% sort){
    # print(year)
    samplesthatyear <- hetdf %>% filter(Year == year)
    
    allelecounts <- samplesthatyear$Allele %>% table %>% data.frame() %>% mutate(Proportion = Freq/sum(Freq))
    # names(allelecounts) <- c("Allele", "Freq", "Proportion")
    names(allelecounts)[1] <- "Allele"
    
    if (length(allelecounts$Allele) < length(allAlleles)){ # One or more of the alleles have frequency 0 that year
      # What alleles are missing?
      missing <- allAlleles[which(!(allAlleles %in% allelecounts$Allele))]
      allelecounts <- rbind(allelecounts, cbind(Allele = missing, Freq = 0, Proportion = 0))
    }    
        
    if (allelecounts$Freq %>% as.numeric() %>% sum() >= minsample){
      PopAlleleCount <- rbind(PopAlleleCount, cbind(year, allelecounts))
    }
  }
  PopAlleleCount <- PopAlleleCount %>% filter(year >= 1900) %>% na.omit()
  PopAlleleCount$Proportion <- as.numeric(PopAlleleCount$Proportion)
  PopAlleleCount$Freq <- as.numeric(PopAlleleCount$Freq)
  
  return(PopAlleleCount)
}

## Function to plot the het genes through time
hetintime <- function(hetcounts, legend = FALSE, title = "Het-z"){
  hetplot <- ggplot(hetcounts, aes(x = year, y = Proportion, colour = Allele)) + 
    xlab("Year") + labs(title = title) +
    # theme(plot.title = element_text(size = rel(1.5), face="bold.italic", family = "Helvetica", hjust = 0.5)) + # I lose control if I use theme_cowplot()
    geom_point(aes(x = year, size = Freq), show.legend = legend) +
    scale_size_continuous(breaks = c(0,5,10,15)) + # Fix the legend circle size
    guides(size=guide_legend(title="Sample\nsize")) +
    geom_line(aes(linetype = Allele)) +
    ylim(0,1) +
    scale_size_area() + # ensures that a value of 0 is mapped to a size of 0
    theme_cowplot() + # simple cowplot theme
    background_grid() # always place this after the theme
  return(hetplot)
}

## Plot the Fst values with its confidence intervals (pre-calculated)
fst.plotter <- function(PopData, FstConfIntv, thistitle = "Chromosome 2", thischr = "chromosome_2", minlen = 5000){
  PopDataFiltered <- PopData %>% filter(variable == "Fst")
  PopDataFiltered$value[which(PopDataFiltered$winlens <= minlen )] <- NA 
  
  p <- ggplot(PopDataFiltered, aes(x = position, y = value, colour = variable)) +
    # geom_ribbon(data = FstConfIntv, aes(ymin = Q1Fst, ymax = maxFst), colour = NA, fill = "darkslategray2", alpha = 0.5) +
    geom_area(data = FstConfIntv, aes(y = maxFst), colour = NA, fill = "darkslategray2", alpha = 0.5) +
    geom_line(alpha = 0.7) +
    theme_bw() +
    labs(title = thistitle) +
    scale_x_continuous(expand = c(0.02, 0.02)) + # Remove the white space on the sides
    theme(plot.title = element_text(size = rel(1.3), face="bold", family = "Helvetica", hjust = 0.5),
          axis.title.y=element_blank(),
          # axis.title.y= element_text(size = rel(1.5)),
          axis.title.x=element_blank(),
          axis.text = element_text(size = rel(1)),
          legend.key=element_blank(), # remove background of the key in the legend
          legend.title=element_blank(),
          legend.text=element_text(size=rel(0.9)),
          legend.position=c(0.9, 0.56),
          legend.spacing.y = unit(0.005, 'cm'), # Put the items of the legend closer to each other
          legend.background=element_blank()) +
    scale_color_manual(values= c("deepskyblue3"), guide= FALSE) + # change color of lines
    ylim(0, 1) +
    geom_line(aes(linetype = variable), alpha = 0.5) + scale_linetype_manual(values=c("solid"), guide = FALSE) +
    geom_point(data = genes %>% filter(chr == thischr) %>% filter(locus != "nwd3"), aes(x = position, shape = locus, fill = locus), colour = "black", size = 3) +
    scale_shape_manual(values= genes[order(genes$locus),] %>% filter(chr == thischr) %>% .$shapes) # Give it specific shapes (but re-order the dataframe so it maches)
  return(p)
}

## Plot the Dxy values
dxy.plotter <- function(PopData, thistitle = "Chromosome 2", thischr = "chromosome_2", legend = "legend", minlen = 5000){
  PopDataFiltered <- PopData %>% filter(variable == "Dxy")
  PopDataFiltered$value[which(PopDataFiltered$winlens <= minlen )] <- NA 
  
  p <- ggplot(PopDataFiltered, aes(x = position, y = value, colour = variable)) +
    geom_area(alpha = 0.5, fill = "dodgerblue4") +
    geom_line(alpha = 0.7) +
    theme_bw() +
    labs(title = thistitle) +
    scale_x_continuous(expand = c(0.02, 0.02)) + # Remove the white space on the sides
    theme(plot.title = element_text(size = rel(1.3), face="bold", family = "Helvetica", hjust = 0.5),
          axis.title.y=element_blank(),
          # axis.title.y= element_text(size = rel(1.5)),
          axis.title.x=element_blank(),
          axis.text = element_text(size = rel(1)),
          legend.key=element_blank(), # remove background of the key in the legend
          legend.title=element_blank(),
          legend.text=element_text(size=rel(0.9)),
          legend.position=c(0.9, 0.56),
          legend.spacing.y = unit(0.005, 'cm'), # Put the items of the legend closer to each other
          legend.background=element_blank()) +
    scale_color_manual(values= c("dodgerblue4"), guide=legend) + # change color of lines
    ylim(0, maxPi) +
    geom_line(aes(linetype = variable), alpha = 0.5) + scale_linetype_manual(values=c("solid"), guide=legend) +
    geom_point(data = genes %>% filter(chr == thischr) %>% filter(locus != "nwd3"), aes(x = position, shape = locus, fill = locus), colour = "black", size = 3) +
    scale_shape_manual(values= genes[order(genes$locus),] %>% filter(chr == thischr) %>% .$shapes) # Give it specific shapes (but re-order the dataframe so it maches)
  return(p)
}

## Plot the Tajima's D values
tajima.plotter <- function(PopData, thistitle = "Chromosome 2", thischr = "chromosome_2", legend = "legend", minlen = 5000){
  PopDataFiltered <- PopData %>% filter(variable == "Tajima" & Pop == "All")
  PopDataFiltered$value[which(PopDataFiltered$winlens <= minlen )] <- NA 
  
  p <- ggplot(PopDataFiltered, aes(x = position, y = value, colour = variable, fill = variable)) + 
    geom_line(alpha = 0.4) + 
    theme_bw() + # Remove background color
    labs(title = thistitle) + 
    scale_x_continuous(expand = c(0.02, 0.02)) + # Remove the white space on the sides
    theme(plot.title = element_text(size = rel(1.5), face="bold", family = "Helvetica", hjust = 0.5), 
          # plot.margin = margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"),
          axis.title.y= element_blank(),
          axis.title.x= element_blank(),
          axis.text = element_text(size = rel(1)),
          legend.key=element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(size=rel(0.9)),
          legend.position=c(0.87, 0.75),  # xy
          legend.background=element_blank()) +
    scale_color_manual(values="darkorange3", labels= "Tajima's D", guide=legend) +
    geom_area(alpha = 0.5, fill = "darkorange3") +
    scale_y_continuous(labels=scaleFUN, limits = c(-3, 5)) +
    geom_point(data = genes %>% filter(chr == thischr) %>% filter(locus != "nwd3"), aes(x = position, shape = locus, fill = locus), colour = "black", size = 3, 
               show.legend = TRUE) + guides(fill=FALSE) + # The guides remove part of the legend
    scale_shape_manual(values= genes[order(genes$locus),] %>% filter(chr == thischr) %>% .$shapes) # Give it specific shapes (but re-order the dataframe so it maches)
  return(p)
}

# Produce a table with the mean values of of a stat
meanstatchr <- function(PopDataChr, stat = "Pi"){
  testgroups <- group_by(PopDataChr %>% filter(variable == stat), Pop)
  pimeans <- summarize(testgroups,
            count = n(),
            meanstat = mean(value, na.rm = TRUE),
            meanwinlens = mean(winlens, na.rm = TRUE))
  return(pimeans)
}

# ============================
### ------Plotting Pi and Theta_W together in one plot -----
# ============================
## Create the plots of each chromosome
PopDataChr1plots <- allpitajima(PopDataChr1, thischr = "chromosome_1", thistitle = "Chromosome 1") # The top in main plot
PopDataChr2plots <- allpitajima(PopDataChr2, thischr = "chromosome_2", thistitle = "Chromosome 2") # The top in supplementary plot
PopDataChr3plots <- allpitajima(PopDataChr3, thischr = "chromosome_3", thistitle = "Chromosome 3", legend = FALSE)
PopDataChr4plots <- allpitajima(PopDataChr4, thischr = "chromosome_4", thistitle = "Chromosome 4", legend = FALSE)
PopDataChr5plots <- allpitajima(PopDataChr5, thischr = "chromosome_5", thistitle = "Chromosome 5", legend = FALSE)
PopDataChr6plots <- allpitajima(PopDataChr6, thischr = "chromosome_6", thistitle = "Chromosome 6", legend = FALSE)
PopDataChr7plots <- allpitajima(PopDataChr7, thischr = "chromosome_7", thistitle = "Chromosome 7", legend = FALSE)

## Process het genes data
hetZcounts <- cbind(yearlycounts(alleledata %>% filter(Het == "Het.Z")), Het = "Het.Z")
hetCcounts <- cbind(yearlycounts(alleledata %>% filter(Het == "Het.c_phen")), Het = "Het.c_phen")
hetScounts <- cbind(yearlycounts(alleledata %>% filter(Het == "Het.s")), Het = "Het.s")

# Rename the S alelles
hetScounts$Allele <- hetScounts$Allele %>% gsub("Big_S", "S", .) %>% gsub("small_s", "s", .) %>% as.factor()

## Plot the het genes historical record
hetz <- hetintime(hetZcounts, legend = TRUE, title = substitute(paste(italic('Het-z'), " (chromosome 1)")))
hetc <- hetintime(hetCcounts, title = substitute(paste(italic('Het-c'), " (chromosome 3)")))
hets <- hetintime(hetScounts, title = substitute(paste(italic('Het-s'), " (chromosome 3)")))

# Put them together in a single object, aligned to each other (important to deal with the legend)
hetplots <- plot_grid(hetz, hetc, hets, ncol = 1, align = "v")

## ----- Just chromosome 1 and 3 in the main figure (MAIN FIGURE) -----
diversitystats2 <- plot_grid(PopDataChr1plots,
                             NULL,
                             PopDataChr3plots,
                             rel_heights = c(1, 0.05, 1),
                             label_size = 12, ncol = 1) 

plot_grid(diversitystats2, NULL, hetplots, ncol = 3, labels = c('A', '', 'B'),
          label_size = 25, rel_widths = c(2, 0.05, 1))
ggsave(snakemake@output$main2chrs, width = 42, height = 28, units = "cm")

plot_grid(PopDataChr2plots,
          NULL,
          PopDataChr4plots,
          NULL,
          PopDataChr5plots,
          NULL,
          PopDataChr6plots,
          NULL,
          PopDataChr7plots,
          rel_heights = c(1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1),
          label_size = 12, ncol = 1) # , labels = c('A', 'B')

ggsave(snakemake@output$mainSupchrs, width = 25, height = 56, units = "cm")

# ============================
### ------ The allele frequencies through time for het-v and het-r -----
# ============================
## Process het genes data
hetVcounts <- cbind(yearlycounts(alleledata %>% filter(Het == "Het.v")), Het = "Het.v")
hetRcounts <- cbind(yearlycounts(alleledata %>% filter(Het == "Het.r_phen")), Het = "Het.r")

# Rename the r alelles
hetRcounts$Allele <- hetRcounts$Allele %>% gsub("Big_R", "R", .) %>% gsub("small_r", "r", .) %>% as.factor()

## --- Plot the genotypes of hetv-hetr through time --- (THE BEST FIGURE IN THE WORLD)
# Prepare a data frame
hetrecordVR <- hetrecord %>% select(-c(Het.Z, Het.c, Het.c_phen, Het.s, Herbivore)) %>% na.omit() # Only het-v and het-r, remove missing data
hetrecordVR$Het.r_phen <- hetrecordVR$Het.r_phen %>% gsub("Big_R", "R", .) %>% gsub("small_r", "r", .) # Replace the het-r alleles for pretty names
hetrecordVR <- hetrecordVR %>% mutate(genotype = paste(Het.r_phen, Het.v, sep='')) %>% select(SampleID, Year, genotype) # Add a column with the genotype, and drop the single genes

# Genotypes alone
genotypes <- melt(hetrecordVR, id = c("SampleID", "Year"))
names(genotypes)[3:4] <- c("Het", "Allele") # to make it compatible with yearlycounts()

genotypesplot <- hetintime(yearlycounts(genotypes), legend = TRUE, title = substitute(paste(italic('het-v'), " and ", italic('het-r'), " genotypes"))) + 
  scale_color_manual(name = 'Genotype', values=c("#B66DE9", "gray", "darkseagreen3", "black")) +
  scale_linetype_manual(name = 'Genotype', values = c(1,3,2,4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_size_area()
  # scale_size(range = c(0, 25))

ggsave(snakemake@output$hetrvtime, width = 5, height = 4)

# -----------
###### Goodness-of-fit
# Are these observation significant?
# The problem with usin Fisher exact test is that it's expecting a 2x2 table, and we know that one genotype is Lethal.
matrixVR<- with(hetrecord %>% select(Het.v,Het.r_phen) %>% na.omit(), table(Het.v, Het.r_phen)) %>% as.matrix
fisher.test(matrixVR, alternative = "two.sided") # p-value < 2.2e-16 # Table produces some weird extra columns with the column names...

matrixVZ<- with(hetrecord %>% select(Het.v, Het.Z) %>% na.omit(), table(Het.v, Het.Z)) %>% as.matrix
fisher.test(matrixVZ, alternative = "two.sided")

matrixVS <- with(hetrecord %>% select(Het.v, Het.s) %>% na.omit(), table(Het.v, Het.s)) %>% as.matrix
fisher.test(matrixVS, alternative = "two.sided")

## But we only have 3 possible genotypes with V-R
chisq.test(matrixVR) # p-value < 2.2e-16

#### Let's use a chi-square but assuming HW as the expected count
# Calculate the frequencies of alleles (counting from the matrix ensures that rows with NAs are not counted)
vr <- c(Vr = matrixVR[3], V1r = matrixVR[4], V1R = matrixVR[2])
freqr <- (vr['Vr'] + vr['V1r'])/ sum(vr)
freqR <- vr['V1R']/ sum(vr)
freqV <- vr['Vr']/ sum(vr)
freqV1 <- (vr['V1R'] + vr['V1r'])/ sum(vr)

# Since VR is letal, that combination is lost
hwp <- c(freqV*freqr, freqV1*freqr, freqV1*freqR)

# *** Is it significant? ***
chisq.test(vr, p = hwp/sum(hwp)) # p-value = 2.336e-11
chisq.test(vr, p = hwp/sum(hwp))$expected #expected values
# With p-value computed with a Monte Carlo test (Hope, 1968) with B replicates
chisq.test(vr, p = hwp/sum(hwp), simulate.p.value=TRUE,B=10000) #  p-value = 9.999e-05

# Make a function for the HW expectations
chisqmyhet <- function(genotypematrix){
  # A is columns, B is rows
  freqA <- sum(genotypematrix[1,])/sum(genotypematrix) #V
  freqa <- sum(genotypematrix[2,])/sum(genotypematrix) #V1
  freqB <- sum(genotypematrix[,1])/sum(genotypematrix) #S
  freqb <- sum(genotypematrix[,2])/sum(genotypematrix) #s
  
  # Expected frequencies under Hardy-Weinberg
  p <- c(freqA*freqB, freqa*freqB, freqA*freqb, freqa*freqb)
  
  return(chisq.test(as.vector(genotypematrix), p = p))
}

chisqmyhet(matrixVZ) # p-value = 0.9236 
chisqmyhet(matrixVS) # p-value = 0.09431

# ============================
### ------ Plotting Fst values -----
# ============================
# Create common xlabel
x.grob <- textGrob("Position (bp)", gp=gpar(fontface="bold", col="black", fontsize=16))
y.grob <- textGrob("Fst", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)

# All chromosomes together
allFstp <- plot_grid(fst.plotter(PopDataChr1, FstConfIntvChr1, thistitle = "Chromosome 1", thischr = "chromosome_1"),
          fst.plotter(PopDataChr2, FstConfIntvChr2, thistitle = "Chromosome 2", thischr = "chromosome_2"),
          fst.plotter(PopDataChr3, FstConfIntvChr3, thistitle = "Chromosome 3", thischr = "chromosome_3"),
          fst.plotter(PopDataChr4, FstConfIntvChr4, thistitle = "Chromosome 4", thischr = "chromosome_4"),
          fst.plotter(PopDataChr5, FstConfIntvChr5, thistitle = "Chromosome 5", thischr = "chromosome_5"),
          fst.plotter(PopDataChr6, FstConfIntvChr6, thistitle = "Chromosome 6", thischr = "chromosome_6"),
          fst.plotter(PopDataChr7, FstConfIntvChr7, thistitle = "Chromosome 7", thischr = "chromosome_7"),
          align = "v", ncol = 1) 
allFst <- grid.arrange(arrangeGrob(allFstp, left = y.grob, bottom = x.grob))

ggsave(plot = allFst, snakemake@output$mainfst, width = 23, height = 34, units = "cm")

# ============================
### ------ Plotting Dxy values -----
# ============================
# Create common xlabel
x.grob <- textGrob("Position (bp)", gp=gpar(fontface="bold", col="black", fontsize=16))
y.grob <- textGrob("Dxy", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)

## Plot the Fst panels
mainDxyp <- plot_grid(
  dxy.plotter(PopDataChr1, thistitle = "Chromosome 1", thischr = "chromosome_1", legend = FALSE),
  dxy.plotter(PopDataChr2, thistitle = "Chromosome 2", thischr = "chromosome_2", legend = FALSE),
  dxy.plotter(PopDataChr3, thistitle = "Chromosome 3", thischr = "chromosome_3", legend = FALSE),
  dxy.plotter(PopDataChr4, thistitle = "Chromosome 4", thischr = "chromosome_4", legend = FALSE),
  dxy.plotter(PopDataChr5, thistitle = "Chromosome 5", thischr = "chromosome_5", legend = FALSE),
  dxy.plotter(PopDataChr6, thistitle = "Chromosome 6", thischr = "chromosome_6", legend = FALSE),
  dxy.plotter(PopDataChr7, thistitle = "Chromosome 7", thischr = "chromosome_7", legend = FALSE),
  align = "v", ncol = 1) 
# Add axes labels to plot
mainDxy <- grid.arrange(arrangeGrob(mainDxyp, left = y.grob, bottom = x.grob))
# Save plot
ggsave(plot = mainDxy, snakemake@output$supdxy, width = 23, height = 34, units = "cm")

# ============================
### ------ Tajima all chromosomes -----
# ============================

# Create common xlabel
x.grob <- textGrob("Position (bp)", gp=gpar(fontface="bold", col="black", fontsize=16))
y.grob <- textGrob("Tajima's D", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)

## Plot the Fst panels
mainTajimap <- plot_grid(
  tajima.plotter(PopDataChr1, thistitle = "Chromosome 1", thischr = "chromosome_1", legend = FALSE),
  tajima.plotter(PopDataChr2, thistitle = "Chromosome 2", thischr = "chromosome_2", legend = FALSE),
  tajima.plotter(PopDataChr3, thistitle = "Chromosome 3", thischr = "chromosome_3", legend = FALSE),
  tajima.plotter(PopDataChr4, thistitle = "Chromosome 4", thischr = "chromosome_4", legend = FALSE),
  tajima.plotter(PopDataChr5, thistitle = "Chromosome 5", thischr = "chromosome_5", legend = FALSE),
  tajima.plotter(PopDataChr6, thistitle = "Chromosome 6", thischr = "chromosome_6", legend = FALSE),
  tajima.plotter(PopDataChr7, thistitle = "Chromosome 7", thischr = "chromosome_7", legend = FALSE),
  align = "v", ncol = 1) 
# Add axes labels to plot
mainTajima <- grid.arrange(arrangeGrob(mainTajimap, left = y.grob, bottom = x.grob))
# Save plot
ggsave(plot = mainTajima, snakemake@output$extratajima, width = 23, height = 45, units = "cm")
# ggsave(plot = mainTajima, snakemake@output[[11]], width = 23, height = 45, units = "cm")
# ggsave(plot = mainTajima, "/Users/Lorena/Dropbox/PhD_UU/Analyses/SnakePipelines/9_DiversityStats/results/figures/Tajima_supp.pdf", width = 23, height = 45, units = "cm")

# ============================
### ------ Plotting Pi per chromosome -----
# ============================

# Make a big table per chromosome for the pi values
pistats <- rbind(cbind(chr = "Chromosome 1", meanstatchr(PopDataChr1, "Pi")),
                 cbind(chr = "Chromosome 2", meanstatchr(PopDataChr2, "Pi")),
                 cbind(chr = "Chromosome 3", meanstatchr(PopDataChr3, "Pi")),
                 cbind(chr = "Chromosome 4", meanstatchr(PopDataChr4, "Pi")),
                 cbind(chr = "Chromosome 5", meanstatchr(PopDataChr5, "Pi")),
                 cbind(chr = "Chromosome 6", meanstatchr(PopDataChr6, "Pi")),
                 cbind(chr = "Chromosome 7", meanstatchr(PopDataChr7, "Pi")) )

# Reverse order of chromosomes, so it's plot nicely in descending order
pistats <- pistats %>% mutate(chr = factor(chr, levels = rev(levels(chr))))

## To get average values
cat(paste0("Average Pi for all samples: ", pistats %>% filter(Pop == "All") %>% .$meanstat %>% mean(),"\n"))
cat(paste0("Average Pi for V samples: ", pistats %>% filter(Pop == "V") %>% .$meanstat %>% mean(),"\n"))
cat(paste0("Average Pi for V1 samples: ", pistats %>% filter(Pop == "V1") %>% .$meanstat %>% mean(),"\n"))

# ============================
### ----- Regions of high Tajima'sD -----
# ============================
# The high Tajima's D can help discover new het genes!
highTajima <- rbind(cbind(chr = "chromosome_1", PopDataChr1 %>% filter(variable == "Tajima" & value >= 2)),
                    cbind(chr = "chromosome_2", PopDataChr2 %>% filter(variable == "Tajima" & value >= 2)),
                    cbind(chr = "chromosome_3", PopDataChr3 %>% filter(variable == "Tajima" & value >= 2)),
                    cbind(chr = "chromosome_4", PopDataChr4 %>% filter(variable == "Tajima" & value >= 2)),
                    cbind(chr = "chromosome_5", PopDataChr5 %>% filter(variable == "Tajima" & value >= 2)),
                    cbind(chr = "chromosome_6", PopDataChr6 %>% filter(variable == "Tajima" & value >= 2)),
                    cbind(chr = "chromosome_6", PopDataChr6 %>% filter(variable == "Tajima" & value >= 2)))

write.csv(highTajima, file=snakemake@output$tajimatable) # Write the table somewhere in case I need it
