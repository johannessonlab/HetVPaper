#!/usr/bin/env Rscript

### PaLDheatmap:Linkage disequilibrium in *Podospora anserina*
#############################################################################
# Part of the Snakemake pipeline PaLD.smk
# Inspired on Hench et al. (2019) Inter-chromosomal coupling between vision and 
# pigmentation genes during genomic divergence, Nature Ecology & evolution 

# The aim of this script is to plot a heatmap of LD from all chromosomes, 
# and between chromosomes 5 and 2
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-09-03
# Version 1
# =======================================

library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(cowplot)
library(wesanderson) # Cute colors

sessionInfo()

# ===========
# Set variables
# ===========
# Pretty colors
pal <- wes_palette("Cavalcanti1", 100, type = "continuous")

## VARIABLES
hetvSTART <- 725000
hetvEND <- 1725000

hetrSTART <- 2250000
hetrEND <- 3250000

buffer <- 100000
sizehetvregion <- 1000000

hetRcoord <- 2683360.5
hetDcoord <- 3056777
hetVcoord <- 1294132.5 - hetvSTART + hetrEND + buffer # Adjusted to the coordinates of the two together
Spok2coord <- 1094363 - hetvSTART + hetrEND + buffer # Adjusted to the coordinates of the two together
Psk7coord <- 900387.5 - hetvSTART + hetrEND + buffer # Adjusted to the coordinates of the two together

MINR2 = 0.000 # chromosomes
MINR2hetrv = 0 # hetrv comparison

midcolor <- 'white'
highcolor <- 'red4'

# ============================
# Prepare functions and variables
# ============================
## Function to prepare the data for plotting
df2ld <- function(chrld){
  # Fix the positions as numbers
  chrld$POS1 <- as.numeric(as.character(chrld$POS1))
  chrld$POS2 <- as.numeric(as.character(chrld$POS2))
  chrld$R.2 <- as.numeric(as.character(chrld$R.2))
  
  # Add a column with the distance between SNPs
  chrld <- chrld %>% mutate(DIST = (POS2 - POS1), MID = POS1 + DIST/2) %>% select(-c(N_CHR, D, Dprime)) # remove columns not-used
  return(chrld)
}

## Function to plot LD triangles of intra- and interchromosomal comparsions
intraldplot <- function(chrld, chrname = "Chromosome 1", filter = MINR2){
  pal <- wes_palette("Zissou1", 50, type = "continuous")
  relevantgenes <- genes %>% filter(CHR == chrld$CHR[1] %>% as.character())
  
  p <- ggplot(chrld %>% filter(R.2 > filter), aes(x = MID, y = DIST)) + geom_point(aes(colour = R.2), shape = 16, alpha = 0.3, size = 0.1) + 
    # scale_colour_gradientn(colours = pal) + # blue 
    scale_color_gradient2(mid = midcolor, high = highcolor) + 
    # theme_cowplot() + # Remove background color, axes, everything
    theme_dark() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.line.y=element_blank(), # remove axis line
          axis.ticks.y=element_blank()) +
    xlim(1, chrdata %>% filter(CHR == chrld$CHR[1]%>% as.character()) %>% .$len) +
    xlab(chrname) +
    geom_point(data = relevantgenes, aes(x = MID, y = 0, shape = Locus, fill = Locus), colour = "black", size = 3) +
    scale_shape_manual(values= relevantgenes[order(relevantgenes$Locus),] %>% .$shapes %>% unique) + # Give it specific shapes (but re-order the dataframe so it maches)
    scale_y_continuous(expand = c(0.02, 0.02)) + 
    # guides(color = FALSE) +
    labs(colour=expression(r^2)) +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
  # theme(legend.position = legendpos)
  return(p)
}

# ============================
# Reading the data
# ============================
cat("*** Reading data... ***\n")
## Data
hetVRld <- read.table(snakemake@input[[1]], header = TRUE, sep = '\t') %>% df2ld()

# Intrachrs
chr1ld.all <- read.table(snakemake@input[[2]], header = TRUE, sep = '\t') %>% df2ld()
chr2ld.all <- read.table(snakemake@input[[3]], header = TRUE, sep = '\t') %>% df2ld()
chr3ld.all <- read.table(snakemake@input[[4]], header = TRUE, sep = '\t') %>% df2ld()
chr4ld.all <- read.table(snakemake@input[[5]], header = TRUE, sep = '\t') %>% df2ld()
chr5ld.all <- read.table(snakemake@input[[6]], header = TRUE, sep = '\t') %>% df2ld()
chr6ld.all <- read.table(snakemake@input[[7]], header = TRUE, sep = '\t') %>% df2ld()
chr7ld.all <- read.table(snakemake@input[[8]], header = TRUE, sep = '\t') %>% df2ld()

## ---- Metadata
## The coordinates of the genes, used by the plotting functions
genesall <- data.frame(Locus = c("Spok2", "MAT", rep("Centromere", 7), 
                                 "het-z", "het-r", "Het-d", "het-s", "het-c", "nwd1", # Deskalov et al 2012: het-s is next to NWD2
                                 "hnwd1", "Het-e", "mod-A", "het-V", "nwd3", "idi-2",
                                 "hnwd3"),
                       MID = c(1094363, 7345878.5, 4479205, 236468.5, 675115.5, 1236388.5, 2062808.5, 238150, 3562141.5,
                               3957511.5, 2683360.5, 3056777, 208504, 495416, 3983565.5,
                               436378.5, 691354.5, 1568320, 1294132.5, 1392335, 1641319,
                               712318.5),
                       CHR = c("chromosome_5", "chromosome_1","chromosome_1", "chromosome_2", "chromosome_3", "chromosome_4", "chromosome_5", "chromosome_6", "chromosome_7",
                               "chromosome_1", "chromosome_2", "chromosome_2", "chromosome_3", "chromosome_3", "chromosome_3",
                               "chromosome_4", "chromosome_4", "chromosome_4", "chromosome_5", "chromosome_5", "chromosome_5",
                               "chromosome_7"),
                       infig = c(1,1,1,1,1,1,1,1,1,  # Use this to filter out loci
                                 1,1,1,1,1,1,
                                 1,1,0,1,1,1,
                                 1),
                       shapes = c(10, 15, rep(19,7), 
                                  17, 24, 4, 5, 6, 8,
                                  12, 18, 13, 25, 2, 11,
                                  9),
                       value = 0)

genes <- genesall %>% filter(infig == 1)

chrdata <- data.frame(CHR = c("chromosome_1","chromosome_2","chromosome_3","chromosome_4","chromosome_5","chromosome_6","chromosome_7"),
                      len = c(8813524, 5165621, 4137471, 3808395, 4734292, 4264132, 4087160)) 

# ============================
# Plot interchromosome comparison
# ============================
cat("*** Plotting interchromosomal comparison of het-v and her-r ... ***\n")
cat(paste0("   minimum r2 used for plotting: ", MINR2, '\n'))

# To plot chromosome lines
hetrline <- data.frame(x1 = hetrSTART, x2 = hetrEND, y1 = -50000, y2 = -50000)
hetvline <- data.frame(x1 = hetrEND + buffer, x2 = hetrEND + buffer + sizehetvregion, y1 = -50000, y2 = -50000) # 

## Same but with artificial background and removing small r2 values
## define the polygons for the background
edge <- 20000 # Extra bp to add a border to the polygons

triangle_r <- data.frame(id = "het-r",
                         x = c(hetrSTART - edge, hetrSTART + (hetrEND - hetrSTART)/2, hetrEND + edge),
                         y = c(0, (hetrEND - hetrSTART) + edge*2, 0))

triangle_v <- data.frame(id = "het-v",
                         x = c(hetrEND + buffer - edge, hetrEND + buffer + (hetvEND - hetvSTART)/2, hetrEND + buffer + (hetvEND - hetvSTART) + edge/2),
                         y = c(0,(hetvEND - hetvSTART) + edge*2, 0) )

rombo <- data.frame(id = "interchr",
                    x = c(hetrEND + buffer/2, hetrSTART + (hetrEND - hetrSTART)/2 + buffer/2 - edge, hetrEND + buffer/2, hetrEND + buffer/2 + sizehetvregion/2 + edge),
                    y = c(buffer - edge*2, (hetrEND - hetrSTART) + buffer, (hetrEND + sizehetvregion - buffer/2)/2 + edge*2, (hetrEND - hetrSTART) + buffer))
# Put together
triangles <- rbind(triangle_r, triangle_v, rombo)

# plot again with background
interchrs <- ggplot(hetVRld %>% filter(R.2 > MINR2hetrv), aes(x = MID, y = DIST)) +
    # theme_dark() +
    scale_color_gradient2(mid = midcolor, high = highcolor) +
    # Chromosome lines
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray25", size = 3, data = hetrline) +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray25", size = 3, data = hetvline) + 
    # Ticks
    geom_segment(aes(x = min(hetVRld$MID), y = -35500, xend = min(hetVRld$MID), yend = -100000), colour = "black", size = 0.7) + # start chr2
    geom_segment(aes(x = hetrEND, y = -35500, xend = hetrEND, yend = -100000), colour = "black", size = 0.7) +# end chr2
    geom_segment(aes(x = hetrEND + buffer, y = -35500, xend = hetrEND + buffer, yend = -100000), colour = "black", size = 0.7) + # start chr5
    geom_segment(aes(x = hetrEND + buffer + sizehetvregion, y = -35500, xend = hetrEND + buffer + sizehetvregion, yend = -100000), colour = "black", size = 0.7) + # end chr5
    # Labels
    annotate(geom="text", x=hetrSTART, y=-150000, label= as.character(hetrSTART), color="black", size = 2.3) + # label start of chr2
    annotate(geom="text", x=hetrEND - 30000, y=-150000, label= as.character(hetrEND), color="black", size = 2.3) + # label end of chr2
    annotate(geom="text", x=hetrEND + buffer, y=-150000, label= as.character(hetvSTART), color="black", size = 2.3) + # label start of chr5
    annotate(geom="text", x=hetrEND + buffer + sizehetvregion, y=-150000, label= as.character(hetvEND), color="black", size = 2.3) + # label end of chr5
    # Chromosome labels
    annotate(geom="text", x=hetrSTART + 500000, y=-200000, label= "Chromosome 2", color="black", size = 5) +
    annotate(geom="text", x=hetrEND + buffer + 500000, y=-200000, label= "Chromosome 5", color="black", size = 5) +
    # Loci
    annotate(geom="text", x=hetRcoord, y=-110000, label= "het-r", color="black", size = 4) +  # label start of chr2
    annotate(geom="text", x=hetDcoord, y=-110000, label= "het-d", color="black", size = 4) +  # label start of chr2
    annotate(geom="text", x=hetVcoord, y=-110000, label= "het-v", color="black", size = 4)  + # label start of chr2
    annotate(geom="text", x=Spok2coord, y=-110000, label= "Spok2", color="black", size = 4)  + # label start of chr2
    # annotate(geom="text", x=Psk7coord, y=-110000, label= "Psk7", color="black", size = 4)  + # label start of chr2
    geom_segment(aes(x = hetRcoord, y = -35500, xend = hetRcoord, yend = -80000), colour = "darkorange", size = 2) + # end chr2
    geom_segment(aes(x = hetDcoord, y = -35500, xend = hetDcoord, yend = -80000), colour = "darkorange3", size = 2) + # end chr2
    geom_segment(aes(x = hetVcoord, y = -35500, xend = hetVcoord, yend = -80000), colour = "darkorange", size = 2) +  # end chr2 
    geom_segment(aes(x = Spok2coord, y = -35500, xend = Spok2coord, yend = -80000), colour = "darkorange3", size = 2) +
    # geom_segment(aes(x = Psk7coord, y = -35500, xend = Psk7coord, yend = -80000), colour = "darkorange3", size = 2) +  
    # Actual points
    geom_point(aes(colour = R.2), shape = 16, alpha = 0.3, size = 0.1) + # points with higher r2 would be a bit bigger
    # geom_point(aes(colour = R.2, size = log1p(R.2)/30), shape = 16, alpha = 0.3) + scale_size_continuous(range = c(0,1), guide = FALSE) + # points with higher r2 would be a bit bigger
    labs(colour=expression(r^2)) + theme(legend.position = "bottom") +
    theme(legend.position = "bottom",
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
    
ggsave(plot = interchrs, snakemake@output[[1]], width = 20, height = 16, units = "cm")

save(interchrs, snakemake@output[[10]]) # Save the object for plotting together with another thing #NEW

# ============================
# Plot intrachromosomal LD
# ============================
cat("*** Plotting intrachromosomal LD (per chr) ... ***\n")

intraldchr1 <- intraldplot(chr1ld.all, chrname = "Chromosome 1")
intraldchr2 <- intraldplot(chr2ld.all, chrname = "Chromosome 2")
intraldchr3 <- intraldplot(chr3ld.all, chrname = "Chromosome 3")
intraldchr4 <- intraldplot(chr4ld.all, chrname = "Chromosome 4")
intraldchr5 <- intraldplot(chr5ld.all, chrname = "Chromosome 5")
intraldchr6 <- intraldplot(chr6ld.all, chrname = "Chromosome 6")
intraldchr7 <- intraldplot(chr7ld.all, chrname = "Chromosome 7")

ggsave(plot = intraldchr1, snakemake@output[[2]], width = 21, height = 14, units = "cm")
ggsave(plot = intraldchr2, snakemake@output[[3]], width = 21, height = 14, units = "cm")
ggsave(plot = intraldchr3, snakemake@output[[4]], width = 21, height = 14, units = "cm")
ggsave(plot = intraldchr4, snakemake@output[[5]], width = 21, height = 14, units = "cm")
ggsave(plot = intraldchr5, snakemake@output[[6]], width = 21, height = 14, units = "cm")
ggsave(plot = intraldchr6, snakemake@output[[7]], width = 21, height = 14, units = "cm")
ggsave(plot = intraldchr7, snakemake@output[[8]], width = 21, height = 14, units = "cm")

cat("*** ... all chromosomes together now ... ***\n")
### Plot all chromosomes together
allchrld <- rbind(chr1ld.all, chr2ld.all, chr3ld.all, chr4ld.all, chr5ld.all, chr6ld.all, chr7ld.all)

allchrsldplot <- ggplot(allchrld %>% filter(R.2 > MINR2), aes(x = MID, y = DIST)) + geom_point(aes(colour = R.2), shape = 16, size = 0.1, alpha = 0.5) +
  scale_color_gradient2(mid = midcolor, high = highcolor) + 
  theme_dark() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_point(data = genes, aes(x = MID, y = 0, shape = Locus, fill = Locus), size = 2) +
  scale_shape_manual(values= genes[order(genes$Locus),] %>% .$shapes %>% unique) + # Give it specific shapes (but re-order the dataframe so it maches)
  facet_grid(CHR ~ ., scales="free_y") +
  labs(colour=expression(r^2))

ggsave(plot = allchrsldplot, snakemake@output[[9]], width = 30, height = 28, units = "cm")

cat("*** DONE! ***\n")
