#conda activate r4
library(ggplot2)
library(RColorBrewer)
library(scales)
library(colorspace)
library(dplyr)

# 1. 102 Malaysia Long-read data [ Keep all records ]
# 2. ho_1926 Short-read data [ Keep selected East and Southeast Asian samples ]
# 3. PanMLS Short-read data [../../reference/PanMLS461/PanMLS461.select.list]
# 4. EGDP Short-read data [ A single level ../../reference/EGDP/sample_info.txt]
# 5. 2018cell Low-coverage short-read data [ Keep all records ]
# 6. Pacific Array data [ Keep all records initially ]
# 7. SG10K Singapore data [ A single level Peranakan_79.Indian_50]



# Get the actual level0 samples
level3_info <- read.table("/path/to/level3.txt", header=TRUE)
pops <- unique(level3_info$POP)

pop_colors <- c(
  'Sea-Bajau' = '#4169E1', # set1: Malaysia
  "Bajau" = '#4197d8',     # 2018Cell
  "Saluan" = '#6D3F95',     # 2018Cell
  "Luzon" = '#bebada', #EGDP
  "Palawan" = '#DA6FFF', #EGDP
  "Cebu" = '#D8BFD8', #EGDP
  "Mindanao" = '#9370DB', #EGDP
  "Mamanwa" = '#DA70D6', #Pacific
  "Manobo" = '#6A5AAA' #Pacific
  )

# Read PCA Results
setwd("/path/to/plink")

pca_data <- read.table("pcs7panel.chrA.level3.geno001.maf005.ld", header=TRUE)
# Name columns (Usually PC1, PC2 Column contains V3 And V4)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4")
pve <- read.table("pve7panel.chrA.level3.geno001.maf005.ld",header=FALSE)
pve
# V1
# 1 0.05604509
# 2 0.02498153
# 3 0.01791483
# 4 0.01516789
# 5 0.01299796
# 6 0.01243624
# 7 0.01106630
# 8 0.01072338
# 9 0.01038914
# 10 0.01027007

final_data <- merge(pca_data, level3_info, by.x = "IID", by.y = "ID", all.x = TRUE)
library(cowplot)

pdf("7panel.chrA.level3.flashpca.pdf", width=10, height=8)
ggplot(final_data, aes(x=PC1, y=PC2, fill=POP)) +
  geom_point(shape = 21, size = 3) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC1(5.60%)") + ylab("PC2(2.50%)") +
  ggtitle("PCA Plot") +
  scale_fill_manual(values = pop_colors) +  #  Apply custom colors
  #coord_fixed(ratio = 1) + # Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )
dev.off()
#+theme(legend.position = "none") # Hide the legend
# ggplot(final_data, aes(x=PC1, y=PC2, fill=Superpop, shape=Superpop)) +
# geom_point(size = 3, aes(color = Group)) +
# xlab("PC1(20.35%)") + ylab("PC2(9.13%)") +
# ggtitle("PCA Plot") +
# scale_fill_manual(values = pop_colors) + # Apply custom colors
# scale_color_manual(values = c("Malaysia" = "black", "Other" = NA)) +
# scale_shape_manual(values = pop_shapes) + # Apply custom markers
# coord_fixed(ratio = 0.715) + # Use a fixed aspect ratio
# theme_minimal() +
# theme(
# plot.title = element_text(hjust = 0.5, size=16, face="bold"),
# legend.title = element_text(size=12),
# legend.text = element_text(size=10)
#)+
# theme(legend.position = "none") # Hide the legend
# dev.off()


# pca_plot <- ggplot(final_data, aes(x=PC1, y=PC2, fill=Superpop)) +
# geom_point(shape = 21, size = 3, aes(color = Group)) +
# #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
# xlab("PC1(20.35%)") + ylab("PC2(9.13%)") +
# ggtitle("PCA Plot") +
# scale_fill_manual(values = pop_colors) + # Apply custom colors
# scale_color_manual(values = c("Malaysia" = "black", "Other" = NA)) +
# #coord_fixed(ratio = 1) + # Use a fixed aspect ratio
# theme_minimal() +
# theme(
# plot.title = element_text(hjust = 0.5, size=16, face="bold"),
# legend.title = element_text(size=12),
# legend.text = element_text(size=10)
#)
# legend_plot <- get_legend(pca_plot)

# # Save the legend
# ggsave("7panel.chrA.level3.v2.pca.legend.pdf", legend_plot, width = 10, height = 6)



