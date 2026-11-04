# Install packages (If not installed)
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("qqman")) install.packages("qqman")

# Load packages
library(tidyverse)  #  Data manipulation
library(qqman)      #  Draw a Manhattan plot
library(ggplot2)
library(dplyr)

setwd("/path/to/Bajau_vs_Malaylike")
CMS <- read.csv("/path/to/merged_all_clean.txt",sep="\t",header=TRUE)
CMS_data <- CMS[,c("CHR","POS","ID","CMS")]

CMS_data$CHR2 <- gsub("chr", "", CMS_data$CHR)
CMS_data$CHR2 <- as.numeric(CMS_data$CHR2)

quantile(CMS_data$CMS, 0.99, na.rm = TRUE)

pdf("CMS_plot.2.pdf",width = 15,height = 8)
# Set colors and threshold lines
gg_manhattan <- ggplot(CMS_data, aes(x = POS, y = CMS, color = factor(CHR2 %% 2))) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("skyblue3", "darkblue")) +
  geom_hline(yintercept = quantile(CMS_data$CMS, 0.99, na.rm = TRUE), linetype = "dashed", color = "red",show.legend = "99%") +
  geom_hline(yintercept = quantile(CMS_data$CMS, 0.999, na.rm = TRUE), linetype = "dashed", color = "darkred",show.legend = "99.9%") +
  facet_grid(. ~ CHR2, scales = "free_x", space = "free_x") +  #  Separate panels by chromosome
  labs(
    title = "CMS Manhattan Plot",
    x = "Chromosomal Position",
    y = "CMS"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    panel.spacing.x = unit(0.1, "lines")
  )
print(gg_manhattan)
dev.off()

