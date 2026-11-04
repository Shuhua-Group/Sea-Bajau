#conda activate r4
library(ggplot2)
library(RColorBrewer)
library(scales)
library(colorspace)
# Group populations by region

# Read all sample metadata
# Five columns, respectively ID POP Superpop Source Language
data_info <- read.table("/path/to/merge.Malaysia_sample.info.txt", header=TRUE,sep="\t")

region_colors <- c(
  'Half_Senoi-Semoq_Beri' = '#FF7F9F', # This project Malaysia
  'Malay' = '#fc8d59', # This project Malaysia
  'Negrito-Bateq' = '#d53e4f', # This project Malaysia
  'Negrito-Mendriq' = '#FF69B4', # This project Malaysia
  'Quarter_Senoi-Semoq_Beri' = '#FFB6C1', # This project Malaysia
  'Senoi-Semai' = '#9370DB', # This project Malaysia
  'Senoi-Semoq_Beri' = '#b8a3e1', # This project Malaysia
  'Senoi-Temiar' = '#9e70a9', # This project Malaysia
  'Brunei' = '#fdae61', # This project Malaysia
  'Bugis' = '#f46d43', # This project Malaysia
  'Dusun' = '#fca875', # This project Malaysia
  'Iban' = '#fee08b', # This project Malaysia
  'Jawa' = '#e58267', # This project Malaysia
  'Kimaragang' = '#ff9e6a', # This project Malaysia
  'Sea-Bajau' = '#4169E1', # This project Malaysia
  'CheWong' = '#F7DD78', # Malaysia
  'Jakun' = '#8FAABB', # Malaysia
  'Mendriq' = '#E8A761', # Malaysia
  'NGS-Dusun' = '#fca800',  # Malaysia
  'Temiar' = '#6a3d9a', # Malaysia
  'Bateq' = '#bebada' # Malaysia
)


pop_shapes <- c(
  'Dusun' = 19, # This project Malaysia
  'Brunei' = 19, # This project Malaysia
  'Bugis' = 19, # This project Malaysia
  'Sea-Bajau' = 19, # This project Malaysia
  'Senoi-Semai' = 19, # This project Malaysia
  'Senoi-Semoq_Beri' = 19, # This project Malaysia
  'Senoi-Temiar' = 19, # This project Malaysia
  'Half_Senoi-Semoq_Beri' = 19, # This project Malaysia
  'Quarter_Senoi-Semoq_Beri' = 19, # This project Malaysia
  'Iban' = 19, # This project Malaysia
  'Jawa' = 19, # This project Malaysia
  'Kimaragang' = 19, # This project Malaysia
  'Malay' = 19, # This project Malaysia
  'Negrito-Bateq' = 19, # This project Malaysia
  'Negrito-Mendriq' = 19, # This project Malaysia
  'Jakun' = 17, # Malaysia
  'Bateq' = 17, # Malaysia
  'CheWong' = 17, # Malaysia
  'Mendriq' = 17, # Malaysia
  'Temiar' = 17, # Malaysia
  'NGS-Dusun' = 17 # Malaysia
)

region_pops <- unique(data_info$Superpop)
subpops <- unique(data_info$POP)

#names(region_colors) <- region_pops

population_colors <- list()
for (regionpop in region_pops) {
  # Get the current regionpop subpopulations
  subpop_in_group <- unique(data_info$POP[data_info$Superpop == regionpop])
  # Check for a single subpop
  if (length(subpop_in_group) == 1) {
    # If there is only one subpop, Use directly regionpop color
    subpop_colors <- region_colors[regionpop]
    names(subpop_colors) <- subpop_in_group
  } else {
    # If there are multiple subpop, Generate distinct colors (Light to dark)
    base_color <- region_colors[regionpop]
    subpop_colors <- colorRampPalette(c(lighten(base_color, 0.4), base_color, darken(base_color, 0.4)))(length(subpop_in_group))
    names(subpop_colors) <- subpop_in_group
  }
  # For the current region_pop assign colors to subpopulations
  population_colors[[regionpop]] <- subpop_colors
}


# Map the generated colors to the data
data_info$Color <- unlist(lapply(data_info$POP, function(sp) {
  regionpop <- data_info$Superpop[data_info$POP == sp][1]
  population_colors[[regionpop]][sp]
}))
head(data_info)

# Check all POP Check that each has a matching shape
unmapped_pops <- setdiff(data_info$POP, names(pop_shapes))
if (length(unmapped_pops) > 0) {
  warning(" Below POP No mapping found for shape: ", paste(unmapped_pops, collapse = ", "))
}

# Mapping shape To POP and create a new column
data_info$Shape <- sapply(data_info$POP, function(pop) {
  if (pop %in% names(pop_shapes)) {
    return(pop_shapes[pop])
  } else {
    return(NA)  #  For unmapped values ,  Set to NA
  }
})
head(data_info)


# Define the marker mapping
shape_mapping <- setNames(data_info$Shape, data_info$POP)

# Define the color mapping
color_mapping <- setNames(data_info$Color, data_info$POP)


# Read PCA Results
pca_data <- read.table("/path/to/malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia.snp.geno005.maf005.ld50k.PCA_out.eigenvec", header=FALSE)
# Name columns (Usually PC1, PC2 Column contains V3 And V4)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
eigenvalues <- read.table("/path/to/malay_ho1926_panmls_SG10K.GRCh38.select1792.Malaysia.snp.geno005.maf005.ld50k.PCA_out.eigenval",header=FALSE)
total_variance <- sum(eigenvalues$V1)
explained_variance <- eigenvalues$V1 / total_variance
explained_variance
# [1] 0.24362481 0.11778541 0.06785235 0.05412587 0.04776476 0.04250785
# [7] 0.03834974 0.03692408 0.03457168 0.03370955 0.03237990 0.03132681
# [13] 0.02924038 0.02841088 0.02797223 0.02755540 0.02732621 0.02696257
# [19] 0.02582658 0.02578294
# Read PCA Data

# Merge PCA data and data_info
merged_data <- merge(pca_data, data_info, by.x = "IID", by.y = "ID", all.x = TRUE)
# Inspect merged results
head(merged_data)


# Plot PCA Plot (PC1 vs PC2)

pdf("102Malay.ho_1926.PanMLS.SG10K.region.Malaysia.geno005.maf005.ld50k.pca.custom.pdf", width=10, height=6)
ggplot(merged_data, aes(x=PC1, y=PC2, color=POP, shape = POP)) +
  geom_point(size=3) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC1 (24.36%)") + ylab("PC2 (11.79%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = color_mapping) +  #  Apply custom colors
  scale_shape_manual(values = shape_mapping) +  #  Apply custom markers
  coord_fixed(ratio = 1) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )
dev.off()

library(cowplot)
pdf("102Malay.ho_1926.PanMLS.SG10K.region.Malaysia.geno005.maf005.ld50k.pca.pdf", width=6, height=6)
ggplot(merged_data, aes(x = PC1, y = PC2, color = POP, shape = POP)) +
  geom_point(size = 3) +
  xlab("PC1 (24.36%)") + ylab("PC2 (11.79%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = color_mapping) +  #  Apply custom colors
  scale_shape_manual(values = shape_mapping) +  #  Apply custom markers
  coord_fixed(ratio = 1) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  theme(legend.position = "none")  #  Hide the legend
dev.off()

# Extract the legend
pca_plot <- ggplot(merged_data, aes(x = PC1, y = PC2, color = POP, shape = POP)) +
  geom_point(size = 3) +
  xlab("PC1 (24.36%)") + ylab("PC2 (11.79%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = color_mapping) +  #  Apply custom colors
  scale_shape_manual(values = shape_mapping) +  #  Apply custom markers
  coord_fixed(ratio = 1) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
legend_plot <- get_legend(pca_plot)

# Save the legend

ggsave("102Malay.ho_1926.PanMLS.SG10K.region.Malaysia.geno005.maf005.ld50k.pca.legend_plot.pdf", legend_plot, width = 10, height = 6)

