#!/usr/bin/env Rscript

### PaLDPlot-nls-grs: A pipeline to estimate linkage disequilibrium in *Podospora anserina*
#############################################################################
# Part of the Snakemake pipeline PaLD.smk
# Based on Hench et al. (2019) Inter-chromosomal coupling between vision and 
# pigmentation genes during genomic divergence, Nature Ecology & evolution 
# S10.Rmd
# The non-linear regression is based on Marroni et al. (2011) Tree Genetics & 
# Genomes 7:1011-1023, with code provided and modified from Tusso et al (2019)
# Mol. Biol. Evol.
# Marroni et al. basically calculate the decay of LD with distance using a 
# nonlinear regression (Hill and Weir 1988) following the equation 1 of 
# Remington et al. (2001). The code uses a non-linear least squares approach 
# (function nls in R) which basically approximate the non-linear function using 
# a linear one and iteratively try to find the best parameter values.

# The aim of this script is to plot the non-linear regression for all samples
# and the two reproductively isolated groups, to see if there are differences. 
# I also can aim to report back some values of rho.
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-08-12 - 2019-09-03
# Version 1
# =======================================

library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(cowplot)
library(hexbin)
library(stringr) # str_extract
library(grid) # for textGrob()
library(gridExtra) # for grid.arrange()

# ============================
# Define some functions
# ============================
## Function to prepare the data for plotting
df2ld <- function(chrld){
  # Fix the positions as numbers
  chrld$POS1 <- as.numeric(as.character(chrld$POS1))
  chrld$POS2 <- as.numeric(as.character(chrld$POS2))
  chrld$R.2 <- as.numeric(as.character(chrld$R.2))
  
  # Add a column with the distance between SNPs
  chrld <- chrld %>% mutate(DIST = (POS2 - POS1)) 
  return(chrld)
}

## *** Make a vector with the predicted values of r2 given the Hill and Weir (1988) 
# and Remington et al. (2001) model ***
# Function based on the code provided by Tusso et al (2019), which in turn came 
# from Knief et al (2017) Molecular Ecology (2017) 26, 1285-1305, who in turn 
# got it from Marroni et al. (2011) :P
# C corresponds to the product between the genetic distance (bp) and the
# population recombination rate (rho) in n number of haplotypes sampled. The
# population recombination rate was calculated as: rho=4*Ne*c, where c is the
# recombination fraction between sites and Ne is the effective population
# size.
HWnonlinear <- function(chrld.df, n, hwst = 0.00001, pop = "All"){
  chr <- chrld.df[1,1]
  cat(paste0("* HWnonlinear of pop ", pop, " for ", chr , " *\n"))
  # Defined starting point for C (the exact value has little effect)
  HW.st <- c(C = hwst) # Starting point for iterations, called like that because Hill and Weir 
  
  # Define variables
  distance <- chrld.df$DIST
  LD.data <- chrld.df$R.2
  
  # Calculate model 
  HW.nonlinear <- nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))), start=HW.st, control=nls.control(maxiter=100))
  tt <- summary(HW.nonlinear)
  
  # This will give you the final parameter for the regression
  new.rho <- tt$parameters[1]        # this is one number
  cat(paste0("    * Calculated rho: ", new.rho, " *\n"))
  
  # Using new.rho, you calculate the predicted value, base on the regression, which will be the values for the line in the plot
  fpoints <-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
  
  newdf <- cbind(chrld.df, hwr2 = fpoints, pop = pop)
  return(newdf)
}

## Plot LD decay just like Marroni et al (2011) (no raw data)
HNnls_plotter <- function(df2ld, legend = 'none', thistitle = "Chromosome 1"){
  h <-  ggplot(df2ld, aes(x= DIST/1000, y = R.2, colour = pop)) + 
    geom_line(aes(y=hwr2), size=1) + # Marroni et al (2011) non-linear regression
    scale_x_continuous(name = "Distance (kb)", expand = c(0,0)) +
    ylim(0, 1) +
    theme_light() +
    labs(title = thistitle) +
    scale_fill_gradient(name = expression(Count~(log[10])), low = rgb(.9,.9,.9), high = "cadetblue4") +
    theme(legend.position = legend, 
          plot.title = element_text(size = rel(1), face="bold", family = "Helvetica", hjust = 0.5),
          axis.title.y= element_blank(), # Remove the axis labels, since I won't be using them
          axis.title.x= element_blank()
    ) + 
    scale_color_manual(values=c("black", "#B66DE9", "darkseagreen3"))
  
  return(h)
}

## Produce a data frame with the bins of r2 (no statistical fitting, but mean of raw data)
ld2bins <- function(chrld){
  ## Categorise distances into intervals of a fixed length (1 kb intervals) and compute mean r2 within each interval
  ## https://www.biostars.org/p/300381/ user rmf style
  chrld$distc <- cut(chrld$DIST, breaks = seq(from = min(chrld$DIST %>% na.omit()) - 1, to = max(chrld$DIST %>% na.omit()) + 1, by = 100)) #
  
  # Then compute mean and/or median r2 within the blocks
  chrldbins <- chrld %>% filter(R.2 != "NaN") %>% group_by(distc) %>% summarise(mean= mean(R.2), median= median(R.2), na.rm = TRUE)
  
  # A helper step to get mid points of our distance intervals for plotting.
  chrldbins <- chrldbins %>% mutate(start= as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                    end= as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                    mid= start + ((end-start)/2) )
  return(chrldbins)
}

## Plot an LD decay curve using two data frames, one coming from HWnonlinear() and another from ld2bins()
# The plot has raw data like Hench et al (2019), it fits the Marroni et al (2011) non-linear regression, 
# and also the average r2 in bins with a smoothing line 
henchHNnls <- function(df2ld.fit, ld2bins, thistitle = "Chromosome 1", legend = 'none'){
  h <-  ggplot(df2ld.fit, aes(x= DIST/1000, y = R.2)) + 
    scale_x_continuous(name = "Distance (kb)", expand = c(0,0)) +
    geom_hex(bins = 100, aes(fill=log10(..count..))) + scale_y_continuous(name = expression(Linkage~(r^2)), expand = c(0.002,0.002)) +
    theme_light() +
    labs(title = thistitle) +
    geom_hex(bins = 100, aes(fill=log10(..count..))) + # fancy hexagons
    scale_fill_gradient(name = expression(Count~(log[10])), low = rgb(.9,.9,.9), high = "cadetblue4") +
    theme(legend.position = legend, 
          plot.title = element_text(size = rel(1), face="bold", family = "Helvetica", hjust = 0.5),
          axis.title.y= element_blank(), # Remove the axis labels, since I won't be using them
          axis.title.x= element_blank()
    ) + 
    geom_smooth(col='red', se=FALSE, size=0.1, method = 'gam', formula = y ~ s(x, bs = "cs")) + # generalized additive mode smoothing # http://statseducation.com/Introduction-to-R/modules/graphics/smoothing/
    geom_line(aes(y=hwr2), size=1, colour = "cadetblue1") + # Marroni et al (2011) non-linear regression
    geom_point(data=ld2bins, aes(x=mid/1000, y=mean), size=0.1, colour="red", alpha = 0.2)  # Mean of bins of data; divide by 1kb so they are in the same scale
  return(h)
}

# ============================
# Reading the data
# ============================

chr1ld.all <- read.table(snakemake@input[[1]], header = TRUE, sep = '\t') %>% df2ld()
chr2ld.all <- read.table(snakemake@input[[2]], header = TRUE, sep = '\t') %>% df2ld()
chr3ld.all <- read.table(snakemake@input[[3]], header = TRUE, sep = '\t') %>% df2ld()
chr4ld.all <- read.table(snakemake@input[[4]], header = TRUE, sep = '\t') %>% df2ld()
chr5ld.all <- read.table(snakemake@input[[5]], header = TRUE, sep = '\t') %>% df2ld()
chr6ld.all <- read.table(snakemake@input[[6]], header = TRUE, sep = '\t') %>% df2ld()
chr7ld.all <- read.table(snakemake@input[[7]], header = TRUE, sep = '\t') %>% df2ld()

chr1ld.grV <- read.table(snakemake@input[[8]], header = TRUE, sep = '\t') %>% df2ld()
chr2ld.grV <- read.table(snakemake@input[[9]], header = TRUE, sep = '\t') %>% df2ld()
chr3ld.grV <- read.table(snakemake@input[[10]], header = TRUE, sep = '\t') %>% df2ld()
chr4ld.grV <- read.table(snakemake@input[[11]], header = TRUE, sep = '\t') %>% df2ld()
chr5ld.grV <- read.table(snakemake@input[[12]], header = TRUE, sep = '\t') %>% df2ld()
chr6ld.grV <- read.table(snakemake@input[[13]], header = TRUE, sep = '\t') %>% df2ld()
chr7ld.grV <- read.table(snakemake@input[[14]], header = TRUE, sep = '\t') %>% df2ld()

chr1ld.grA <- read.table(snakemake@input[[15]], header = TRUE, sep = '\t') %>% df2ld()
chr2ld.grA <- read.table(snakemake@input[[16]], header = TRUE, sep = '\t') %>% df2ld()
chr3ld.grA <- read.table(snakemake@input[[17]], header = TRUE, sep = '\t') %>% df2ld()
chr4ld.grA <- read.table(snakemake@input[[18]], header = TRUE, sep = '\t') %>% df2ld()
chr5ld.grA <- read.table(snakemake@input[[19]], header = TRUE, sep = '\t') %>% df2ld()
chr6ld.grA <- read.table(snakemake@input[[20]], header = TRUE, sep = '\t') %>% df2ld()
chr7ld.grA <- read.table(snakemake@input[[21]], header = TRUE, sep = '\t') %>% df2ld()

# SAMPLE sizes
n.all <- 106
n.grV <- 57
n.grA <- 49

# ============================
# Analysis
# ============================
cat("*** Fitting the data to the Remington et al. (2001) equation ***\n")
chr1ld.all.fit <- HWnonlinear(chr1ld.all, n = n.all, pop = "All")
chr2ld.all.fit <- HWnonlinear(chr2ld.all, n = n.all, pop = "All")
chr3ld.all.fit <- HWnonlinear(chr3ld.all, n = n.all, pop = "All")
chr4ld.all.fit <- HWnonlinear(chr4ld.all, n = n.all, pop = "All")
chr5ld.all.fit <- HWnonlinear(chr5ld.all, n = n.all, pop = "All")
chr6ld.all.fit <- HWnonlinear(chr6ld.all, n = n.all, pop = "All")
chr7ld.all.fit <- HWnonlinear(chr7ld.all, n = n.all, pop = "All")

# For RI group with the V allele
chr1ld.grV.fit <- HWnonlinear(chr1ld.grV, n = n.grV, pop = "V")
chr2ld.grV.fit <- HWnonlinear(chr2ld.grV, n = n.grV, pop = "V")
chr3ld.grV.fit <- HWnonlinear(chr3ld.grV, n = n.grV, pop = "V")
chr4ld.grV.fit <- HWnonlinear(chr4ld.grV, n = n.grV, pop = "V") 
chr5ld.grV.fit <- HWnonlinear(chr5ld.grV, n = n.grV, pop = "V")
chr6ld.grV.fit <- HWnonlinear(chr6ld.grV, n = n.grV, pop = "V")
chr7ld.grV.fit <- HWnonlinear(chr7ld.grV, n = n.grV, pop = "V")

# For RI group with the V1 allele
chr1ld.grA.fit <- HWnonlinear(chr1ld.grA, n = n.grA, pop = "V1")
chr2ld.grA.fit <- HWnonlinear(chr2ld.grA, n = n.grA, pop = "V1")
chr3ld.grA.fit <- HWnonlinear(chr3ld.grA, n = n.grA, pop = "V1")
chr4ld.grA.fit <- HWnonlinear(chr4ld.grA, n = n.grA, pop = "V1")
chr5ld.grA.fit <- HWnonlinear(chr5ld.grA, n = n.grA, pop = "V1")
chr6ld.grA.fit <- HWnonlinear(chr6ld.grA, n = n.grA, pop = "V1")
chr7ld.grA.fit <- HWnonlinear(chr7ld.grA, n = n.grA, pop = "V1")

# Chromosome 4 is acting weird, likely because of the high Tajima's D blocks at the extremes
chr4ld.all.fit_reduce <- HWnonlinear(chr4ld.all %>% filter(POS1 > 1000000  & POS1 < 3500000), n = n.all, pop = "All")

cat("*** Plotting full LD decay curves for all samples ***\n")
# Extract the legend alone
legend <- cowplot::get_legend(henchHNnls(chr7ld.all.fit, ld2bins(chr7ld.all.fit), "Chromosome 7", legend = 'bottom'))

# Create common x and y labels
x.grob <- textGrob("Distance (Kbp)", gp=gpar(col="black", fontsize=11))
y.grob <- textGrob(expression(Linkage~(r^2)), gp=gpar(fontface="bold", col="black", fontsize=11), rot=90)

plotsLD <- plot_grid(
  henchHNnls(chr1ld.all.fit, ld2bins(chr1ld.all.fit), thistitle = "Chromosome 1"),
  henchHNnls(chr2ld.all.fit, ld2bins(chr2ld.all.fit), thistitle = "Chromosome 2"),
  henchHNnls(chr3ld.all.fit, ld2bins(chr3ld.all.fit), thistitle = "Chromosome 3"),
  henchHNnls(chr4ld.all.fit, ld2bins(chr4ld.all.fit), thistitle = "Chromosome 4"),
  henchHNnls(chr5ld.all.fit, ld2bins(chr5ld.all.fit), thistitle = "Chromosome 5"),
  henchHNnls(chr6ld.all.fit, ld2bins(chr6ld.all.fit), thistitle = "Chromosome 6"),
  henchHNnls(chr7ld.all.fit, ld2bins(chr7ld.all.fit), thistitle = "Chromosome 7"),
  henchHNnls(chr4ld.all.fit_reduce, ld2bins(chr4ld.all.fit_reduce), thistitle = "Chromosome 4:1-3.5Mb"),
  align = "h", ncol = 4, nrow = 2)

# Add axes labels to plot
decayallfits <- grid.arrange(arrangeGrob(plotsLD, left = y.grob, bottom = x.grob))
# Add legend at the bottom
decayallfits.final <- plot_grid(decayallfits, legend, ncol = 1, nrow = 2, rel_heights = c(4, 0.3))

# Save
ggsave(plot = decayallfits.final, snakemake@output[[1]], width = 22, height = 13, units = "cm")

### -------
cat("*** Plotting Remington's curves per reproductively isolated group ***\n")
# Extract the legend alone
legend <- cowplot::get_legend(HNnls_plotter(rbind(chr7ld.all.fit, chr7ld.grV.fit, chr7ld.grA.fit), "Chromosome 7", legend = 'right'))

# Create common x and y labels
x.grob <- textGrob("Distance (Kbp)", gp=gpar(col="black", fontsize=11))
y.grob <- textGrob(expression(Linkage~(r^2)), gp=gpar(fontface="bold", col="black", fontsize=11), rot=90)

decayplots <- plot_grid(
  HNnls_plotter(rbind(chr1ld.all.fit, chr1ld.grV.fit, chr1ld.grA.fit), thistitle = "Chromosome 1"),
  HNnls_plotter(rbind(chr2ld.all.fit, chr2ld.grV.fit, chr2ld.grA.fit), thistitle = "Chromosome 2"),
  HNnls_plotter(rbind(chr3ld.all.fit, chr3ld.grV.fit, chr3ld.grA.fit), thistitle = "Chromosome 3"),
  HNnls_plotter(rbind(chr4ld.all.fit, chr4ld.grV.fit, chr4ld.grA.fit), thistitle = "Chromosome 4"),
  HNnls_plotter(rbind(chr5ld.all.fit, chr5ld.grV.fit, chr5ld.grA.fit), thistitle = "Chromosome 5"),
  HNnls_plotter(rbind(chr6ld.all.fit, chr6ld.grV.fit, chr6ld.grA.fit), thistitle = "Chromosome 6"),
  HNnls_plotter(rbind(chr7ld.all.fit, chr7ld.grV.fit, chr7ld.grA.fit), thistitle = "Chromosome 7"),
  legend,
  align = "h", ncol = 4, nrow = 2) 

# Add axes labels to plot
decayplots.final <- grid.arrange(arrangeGrob(decayplots, left = y.grob, bottom = x.grob))

# Save
ggsave(plot = decayplots.final, snakemake@output[[2]], width = 20, height = 11, units = "cm")

# ============================
# How fast is the LD decaying per chromosome?
# ============================
# If we use a threshold of LD (say, 0.2), what is the maximum average linkage at that point?
midr2threshold <- function(chrld.all.fit, threshold = 0.2){
  ld2binsch <- ld2bins(chrld.all.fit)
  tablethres <- ld2binsch[order(ld2binsch$mean),] %>% filter(mean <= threshold) 
  thresdist <- tablethres$mid %>% min(na.rm = TRUE) # Get the smallest distance value of the ranges lower than the threshold
  thischr <- chrld.all.fit$CHR[1] # assume only one chr present
  
  return(data.frame(chr = thischr, dist = thresdist, threshold = threshold))
}

midr2s <- rbind(midr2threshold(chr1ld.all.fit),
                midr2threshold(chr2ld.all.fit), 
                midr2threshold(chr3ld.all.fit),
                midr2threshold(chr4ld.all.fit),
                midr2threshold(chr5ld.all.fit),
                midr2threshold(chr6ld.all.fit),
                midr2threshold(chr7ld.all.fit) )

# Print it to the log
print(midr2s)
# Save the file
write.csv(midr2s, file=snakemake@output[[3]]) # Write the table somewhere in case I need it
