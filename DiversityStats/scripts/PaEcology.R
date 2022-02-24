#!/usr/bin/env Rscript

### PaEcology: Explore the ecology of the mating groups in P. anserina
#############################################################################
# Part of the Snakemake pipeline DiversityStats.smk
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019-09-12
# Version 1
# =======================================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE)
library(cowplot) 
library(reshape2) # for melt

# ============================
# Reading the data
# ============================

hetrecord <- read.csv("../data/PaAllelesPaper.csv", header = TRUE, na.strings = "NA", sep=",") # Replace the path, as the relative paths won't work
sampling2017 <- read.csv("../data/Sampling2017_20220517.csv", head = TRUE, na.strings = "NA", sep=",") # Replace the path, as the relative paths won't work
sampling2017 <- sampling2017 %>% filter(Indiv_count != "-" & Mating_phenotype != "P. comata") %>% select(c(Wa_numbers, Mating_phenotype, Dung, Indiv_count))
sampling2017 <- sampling2017[!(is.na(sampling2017$Mating_phenotype)), ]  # Remove columns without phenotype data


# ============================
# Are there differences in substrate?
# ============================
countherbis <- with(hetrecord, table(Het.v, Herbivore))
herbiabundance <- countherbis %>% melt

cat(paste("We have substrate data for a total of", herbiabundance$value %>% sum, "strains"))

cat("Is this significant?")
cat("(We ignored sheep due to low sample sizes)")
chisq.test(countherbis[,1:2]) # Not significant p-value = 0.4959
# There is no pattern of a particular het-v allele associating with a specific herbivore (Horse or Rabbit)

# What about individually? null hypothesis that they are equally distributed
chisq.test(countherbis[,1]) # p-value = 0.8886 for Horse
chisq.test(countherbis[,2]) # p-value = 0.2087 for Rabbit

# Plot
substrateplot <- ggplot(herbiabundance, aes(fill=Het.v, x=Herbivore, y=value)) + 
  geom_bar(position="stack", stat="identity") +
  theme_cowplot() + 
  ylab("Count") + ggtitle("1991-2016") +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values= c("#B66DE9", "darkseagreen3"), guide = "none") +
  annotate("text", x = 1, y=52, label = "n.s.") + 
  annotate("text", x = 2, y=32, label = "n.s.") 

ggsave(plot = substrateplot, "../results/figures/Fig6B.pdf", width = 2.5, height = 4)

# ============================
# Do they co-exist in the same dung piece?
# ============================
dungs <- with(sampling2017 , table(Dung, Mating_phenotype)) %>% melt(id = c("Dung"))
names(dungs) <- c("Dung", "Group", "Count")

dungplot <- ggplot(dungs, aes(fill=Group, x=as.factor(Dung), y=Count)) + 
  geom_bar(position="stack", stat="identity") +
  theme_cowplot() +
  xlab("Dung sample") + ggtitle("2017") + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = c(.9, .65)) +
  scale_fill_manual(values= c("#B66DE9","darkseagreen3", "gray"), name = "Inferred\n Group", 
                    labels= c("rV", "RV1"))

ggsave(plot = dungplot, "./results/figures/Fig4C.pdf", width = 8.2, height = 3)
