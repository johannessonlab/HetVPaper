#!/usr/bin/env Rscript

### PaMutantsFertility: Effects of mutants of Pa_5_12720 on fertility
#############################################################################
# Part of the Snakemake pipeline DiversityStats.smk
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2021-12-12
# Version 1
# =======================================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE)
library(cowplot)
library(rstatix)
library(ggpubr)

# ============================
# Reading the data
# ============================
fertilitydata <- read.csv("../data/Mutants_fertility_data.csv", header = TRUE, na.strings = "NA", sep=",") # Replace the path, as the relative paths won't work

# ============================
# Are they significantly different?
# ============================
# https://statsandr.com/blog/anova-in-r/
# The data might not be normally distributed, we don't have large sample sizes, 
# and the groups have different variances.
# library("car")
# leveneTest(Perithecia ~ Focal_strain, data = fertilitydata %>% filter(Condition == "Male_fertility_to_V1R"))

# res_aov <- aov(Perithecia ~ Focal_strain,
#                data = fertilitydata %>% filter(Condition == "Female_fertility_to_V1r"))
# hist(res_aov$residuals)
# library(car)
# car::qqPlot(res_aov$residuals, id = FALSE) # id = FALSE to remove point identification
# shapiro.test(res_aov$residuals)
# we reject the hypothesis that residuals follow a normal distribution (p-value = 0.002228).

# So I'll use a non-parametric test

## --- Female fertility ---
datafemale <- fertilitydata %>% filter(Condition == "Female_fertility_to_V1r")
k1 <- kruskal.test(Perithecia ~ Focal_strain, data = datafemale)
# Kruskal-Wallis chi-squared = 26.383, df = 6, p-value = 0.0001888
# we reject the hypothesis that all means are equal

pairwise.wilcox.test(datafemale$Perithecia, datafemale$Focal_strain, p.adjust.method = 'bonferroni')
# pairwise.wilcox.test(datafemale$Perithecia, datafemale$Focal_strain, p.adjust.method = 'BH') # similar results

# Same as
rstatix::pairwise_wilcox_test(datafemale, Perithecia ~ Focal_strain, p.adjust.method = "bonferroni")

## --- Male fertility ---
datamale <- fertilitydata %>% filter(Condition == "Male_fertility_to_V1R")
k2 <- kruskal.test(Perithecia ~ Focal_strain, data = datamale)
# Kruskal-Wallis chi-squared = 40.862, df = 6, p-value = 3.083e-07
# we reject the hypothesis that all means are equal

pairwise.wilcox.test(datamale$Perithecia, datamale$Focal_strain, p.adjust.method = 'bonferroni')
# pairwise.wilcox.test(datafemale$Perithecia, datafemale$Focal_strain, p.adjust.method = 'BH') # similar results

# ============================
# Plot
# ============================
## Reorder the levels
# VERY important to reorder because the assignments of add_xy_position() follows 
# the native order of levels for Focal_strain of R/ggplot, but ggboxplot follows 
# the order of the original database, so the labels of stats get assigned wrong!
fertilitydata$Focal_strain <- factor(fertilitydata$Focal_strain, 
                                     levels = c("V1r", "Pa_5_12720", "Pa_5_12720 Y233F", "Pa_5_12720 Y233A", "Pa_5_12720+Pa_5_12710", "Pa_5_12720 Y233F+Pa_5_12710","Pa_5_12720 Y233A+Pa_5_12710"))

# Change name of treatments to have white spaces
# We tend to write Vr, V1r, etc, but in the paper we put her-r before het-v...
fertilitydata$Condition <- as.factor(fertilitydata$Condition)
fertilitydata$Condition<- recode_factor(fertilitydata$Condition, Female_fertility_to_V1r = "Female fertility to rV1", Male_fertility_to_V1R = "Male fertility to RV1")

# Calculate the stats
stat.test <- fertilitydata %>%
  group_by(Condition) %>%
  pairwise_wilcox_test(Perithecia ~ Focal_strain, p.adjust.method = "bonferroni")

# Add adjusted p-values
stat.test <- stat.test %>% add_xy_position(x = "Focal_strain")

# Correct coordinates after removing the not significant comparisons
# Otherwise a lot of white space is left
stat.test2 <- stat.test %>% filter(p.adj.signif != "ns")

stat.test.f <- stat.test2 %>% filter(Condition == "Female fertility to rV1")
maxcoord <- max(datafemale$Perithecia)
stat.test.f$y.position <- maxcoord + seq(1, nrow(stat.test.f), 1 )*2

stat.test.m <- stat.test2 %>% filter(Condition == "Male fertility to RV1")
maxcoord <- max(datamale$Perithecia)
stat.test.m$y.position <- maxcoord + seq(1, nrow(stat.test.m), 1 )*2

stat.test2 <- rbind(stat.test.f, stat.test.m)

## Plot with ggpubr
# https://www.datanovia.com/en/blog/ggpubr-how-to-add-adjusted-p-values-to-a-multi-panel-ggplot/
ggboxplot(fertilitydata, x = "Focal_strain", y = "Perithecia", 
  fill = "Focal_strain") +
  geom_jitter(alpha = 0.3, size = 1, width = 0.1) +
  xlab(NULL) +
  ylab("Perithecia per contact") +
  facet_wrap(Condition ~ ., ncol = 1, strip.position="top") +
  rotate_x_text(angle = 60) +
  theme(strip.background = element_rect(fill="gray94"),
        axis.text.x = element_text(size=10), #angle = 60, hjust=1, 
        strip.text.x = element_text(size = 12),
        axis.title.y=element_text(size=14),
        legend.position = "none") +
  scale_fill_brewer(palette="RdYlBu") +
  # scale_x_discrete(labels= c("rV1", "Pa_5_12720", "Pa_5_12720 Y233F", "Pa_5_12720 Y233A", "Pa_5_12720 +\nPa_5_12710", "Pa_5_12720 Y233F +\nPa_5_12710","Pa_5_12720 Y233A +\nPa_5_12710")) +
  scale_x_discrete(labels= c("rV1", "het-Va", "het-Va Y233F", "het-Va Y233A", "het-Va +\nhet-Vb", "het-Va Y233F +\nhet-Vb","het-Va Y233A +\nhet-Vb")) +
  stat_pvalue_manual(stat.test2, label = "p.adj.signif", hide.ns = TRUE, tip.length = 0) 

ggsave("../results/FigS11b_fertility.png", width = 3.3, height = 6.5) # Replace the path, as the relative paths won't work
