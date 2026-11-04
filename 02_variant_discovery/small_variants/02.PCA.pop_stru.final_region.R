#conda activate r4
library(ggplot2)
library(RColorBrewer)
library(scales)
library(colorspace)
# Group populations by region

# Read all sample metadata
# Five columns, respectively ID POP Superpop Source Language
data_info <- read.table("/path/to/merge.sample.info.txt", header=TRUE,sep="\t")

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
  'Bateq' = '#bebada', # Malaysia
  'America' = "#B99B60",
  "EastAsia" = "#62C15D",
  'Indian' = "#D18C7F",
  'Kusunda' = '#FFD700',
  'Kyrgyz' = "#1F77B4",
  'Oceania' = "#1BAC7C",
  'Singapore' = '#00CED1',
  'SouthAsia' = "#E78AC3"
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
  'NGS-Dusun' = 17, # Malaysia
  'Achang' = 0,
  'Ami' = 1,
  'Atayal' = 2,
  'Bai' = 3,
  'Balochi' = 4,
  'Bengali' = 6,
  'Birhor' = 7,
  'Blang' = 8,
  'Blang-DY' = 9,
  'Bougainville' = 10,
  'Bouyei' = 11,
  'Brahui' = 12,
  'Burmese' = 15,
  'Burusho' = 16,
  'Cambodian' = 17,
  'Chaoxian' = 18,
  'Dai' = 17,
  'Daur' = 21,
  'Deang' = 22,
  'Dong' = 23,
  'Dongxiang' = 24,
  'Dongxiang-YZ' = 25,
  'Drung' = 0,
  'Even' = 2,
  'Evenki' = 3,
  'Gelao' = 4,
  'Gin' = 5,
  'Gujarati' = 6,
  'Han' = 8,
  'Han_N' = 9,
  'Hani' = 10,
  'Hawaiian' = 11,
  'Hazara' = 12,
  'Hezhen' = 13,
  'Hui' = 14,
  'Indian' = 16,
  'Japanese' = 18,
  'Jingpo' = 15,
  'Jino' = 21,
  'Kalash' = 22,
  'Kharia' = 23,
  'Kinh' = 25,
  'Korean' = 0,
  'Kusunda' = 1,
  'Kyrgyz' = 2,
  'Lahu' = 3,
  'Lisu' = 4,
  'Lodhi' = 5,
  'Makrani' = 6,
  'Mala' = 7,
  'Man' = 9,
  'Miao' = 11,
  'Mongol' = 12,
  'Mongol-DMQ' = 13,
  'Mongol-TH' = 14,
  'Mongola' = 15,
  'Mosuo' = 16,
  'Mulam' = 10,
  'NAN_Melanesian' = 18,
  'Naxi' = 19,
  'Nu' = 22,
  'Onge' = 23,
  'Oroqen' = 24,
  'POP' = 25,
  'Papuan' = 0,
  'Pathan' = 1,
  'Peranakan_Indian' = 2,
  'Peranakan_Caucasian' = 3,
  'Peranakan_Chinese' = 4,
  'Peranakan_Chinese_Peranakan_Eurasian' = 5,
  'Peranakan_Eurasian' = 6,
  'Pumi' = 7,
  'Punjabi' = 8,
  'Qiang' = 9,
  'Salar' = 11,
  'She' = 16,
  'Sherpa' = 14,
  'Shui' = 18,
  'Sindhi' = 15,
  'Thai' = 21,
  'Tibetan' = 22,
  'Tu' = 23,
  'Tu-HZ' = 24,
  'Tu-XN' = 25,
  'Tujia' = 0,
  'Uygur' = 1,
  'Wa' = 2,
  'Wa-XM' = 3,
  'Xibo' = 4,
  'Yakut' = 5,
  'Yao' = 6,
  'Yi' = 7,
  'Yugur' = 8,
  'Zhuang' = 9
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
pca_data <- read.table("/path/to/malay_ho1926_panmls_SG10K.GRCh38.select1792.snp.geno005.maf005.ld50k.PCA_out.eigenvec", header=FALSE)
# Name columns (Usually PC1, PC2 Column contains V3 And V4)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
eigenvalues <- read.table("/path/to/malay_ho1926_panmls_SG10K.GRCh38.select1792.snp.geno005.maf005.ld50k.PCA_out.eigenval",header=FALSE)
total_variance <- sum(eigenvalues$V1)
explained_variance <- eigenvalues$V1 / total_variance
explained_variance
# [1] 0.43324497 0.11145703 0.07366578 0.06283489 0.04448425 0.03477528
# [7] 0.02837717 0.02318596 0.02077724 0.01992971 0.01783794 0.01756597
# [13] 0.01664502 0.01541173 0.01438353 0.01404105 0.01331397 0.01299228
# [19] 0.01287140 0.01220481
# Read PCA Data

# Merge PCA data and data_info
merged_data <- merge(pca_data, data_info, by.x = "IID", by.y = "ID", all.x = TRUE)
# Inspect merged results
head(merged_data)


# Plot PCA Plot (PC1 vs PC2)

pdf("102Malay.ho_1926.PanMLS.SG10K.region.geno005.maf005.ld50k.pca.custom.pdf", width=20, height=8)
ggplot(merged_data, aes(x=PC2, y=PC1, color=POP, shape = POP)) +
  geom_point(size=3) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC2 (11.46%)") + ylab("PC1 (43.32%)") +
  # xlab("PC1 (43.32%)") + ylab("PC2 (11.46%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = color_mapping) +  #  Apply custom colors
  scale_shape_manual(values = shape_mapping) +  #  Apply custom markers
  coord_fixed(ratio = 1.8) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )
dev.off()

library(cowplot)
pdf("102Malay.ho_1926.PanMLS.SG10K.region.geno005.maf005.ld50k.pca.pdf", width=8, height=8)
ggplot(merged_data, aes(x = PC1, y = PC2, color = POP, shape = POP)) +
  geom_point(size = 3) +
  xlab("PC1 (43.32%)") + ylab("PC2 (11.46%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = color_mapping) +  #  Apply custom colors
  scale_shape_manual(values = shape_mapping) +  #  Apply custom markers
  coord_fixed(ratio = 0.5) +  #  Use a fixed aspect ratio
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
  xlab("PC1 (43.32%)") + ylab("PC2 (11.46%)") +
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

ggsave("102Malay.ho_1926.PanMLS.SG10K.region.geno005.maf005.ld50k.pca.legend_plot.pdf", legend_plot, width = 10, height = 6)

