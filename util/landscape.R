#!/usr/bin/env Rscript
rm(list=ls())

library(ggplot2)
library(tidyverse)
library(RColorBrewer)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: script.R <input_csv> <output_image>")
}

input_tsv <- args[1]
output_image <- args[2]

d <- read.csv(input_tsv, sep = "\t")

dnew <- d[d$Kimura < 60,]

# Define the groups

dnew$SuperFamily <- factor(dnew$SuperFamily, levels=c(
  "ClassI/Unclassified", "DIRS/DIRS", "LINE", "LINE/CR1", "LINE/I", "LINE/L2", "LINE/R1", "LINE/R2",
  "LINE/RTE", "LTR", "LTR/Copia", "LTR/Gypsy", "LTR/LARD", "LTR/TRIM", "SINE","ClassII/Unclassified", "DNA", "DNA/Crypton",
  "DNA/hAT", "DNA/Maverick", "DNA/Merlin", "DNA/MITE", "DNA/P", "DNA/PIF-Harbinger",
  "DNA/PiggyBac", "RC/Helitron"
))

red_gradient <- colorRampPalette(c("firebrick","gold"))(15)
blue_gradient <- colorRampPalette(c("dodgerblue3", "forestgreen"))(10)
green_gradient <- colorRampPalette(c("lightseagreen", "lightseagreen"))(1)
class2grad <- colorRampPalette(c("blue", "blue"))(1)
class1grad <- colorRampPalette(c("red", "red"))(1)

color_mapping <- c(
  "ClassII/Unclassified" = class2grad[1],
  "ClassI/Unclassified" = class1grad[1],
  "DIRS/DIRS" = red_gradient[2],
  "DNA" = blue_gradient[2],
  "DNA/Crypton" = blue_gradient[3],
  "DNA/hAT" = blue_gradient[4],
  "DNA/Maverick" = blue_gradient[5],
  "DNA/Merlin" = blue_gradient[6],
  "DNA/MITE" = blue_gradient[7],
  "DNA/P" = blue_gradient[8],
  "DNA/PIF-Harbinger" = blue_gradient[9],
  "DNA/PiggyBac" = blue_gradient[10],
  "LINE" = red_gradient[3],
  "LINE/CR1" = red_gradient[4],
  "LINE/I" = red_gradient[5],
  "LINE/L2" = red_gradient[6],
  "LINE/R1" = red_gradient[7],
  "LINE/R2" = red_gradient[8],
  "LINE/RTE" = red_gradient[9],
  "LTR" = red_gradient[10],
  "LTR/Copia" = red_gradient[11],
  "LTR/Gypsy" = red_gradient[12],
  "LTR/LARD" = red_gradient[13],
  "LTR/TRIM" = red_gradient[14],
  "RC/Helitron" = green_gradient[1],
  "SINE" = red_gradient[15]
)

p <- ggplot(dnew, aes(x=Kimura, y=Percentage.of.Genome,fill=SuperFamily,group=SuperFamily)) +
  geom_bar(stat="identity",colour="white", size = 0.35) +
  labs(x = "Kimura distance (adjusted CpG model)", y = "% of genome") +
  scale_fill_manual(values = color_mapping) +
  #theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.title = element_text(size = 12),  # Smaller legend title
    legend.text = element_text(size = 10),   # Smaller legend text
    legend.key.size = unit(0.8, "lines")    # Smaller legend key size
  )

ggsave(filename = output_image, plot = p, width = 11, height = 7)
