#conda activate r4
library(ggplot2)
library(RColorBrewer)
library(scales)
library(colorspace)
# Group populations by region

# Read all sample metadata
# Five columns, respectively ID POP Superpop Source Language
data_info <- read.table("/path/to/merge.Malaysia_sample.90_2018Cell.EGDP.EastSouth.info.txt", header=TRUE,sep="\t")

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
  'Bajau' = '#4197d8', # 2018Cell
  'Saluan' = '#6D3F95',  # 2018Cell
  'Ami' = '#8eca7c', # EastSouthAsia
  'Atayal' = '#8cca7c', # EastSouthAsia
  'Birhor' = '#efaacc', # EastSouthAsia
  'Blang' = '#8ac979', # EastSouthAsia
  'Blang-DY' = '#89c979', # EastSouthAsia
  'Bougainville' = '#6dc49e', # EastSouthAsia
  'Cambodian' = '#86c877', # EastSouthAsia
  'Dai' = '#83c775', # EastSouthAsia
  'Kharia' = '#df85b7', # EastSouthAsia
  'NAN_Melanesian' = '#1aac7c', # EastSouthAsia
  'Peranakan_Caucasian' = '#02787a', # EastSouthAsia
  'Peranakan_Chinese' = '#74cdd9', # EastSouthAsia
  'Peranakan_Chinese_Peranakan_Eurasian' = '#03a2a4', # EastSouthAsia
  'Peranakan_Eurasian' = '#57c5cf', # EastSouthAsia
  'Peranakan_Indian' = '#3ec1c8', # EastSouthAsia
  'Thai' = '#449544', # EastSouthAsia
  'Bangladesh' = '#d6a5ab', # EGDP
  'Borneo' = '#D0AABE', # EGDP
  'Cebu' = '#D3C7B9', # EGDP
  'India' = '#7FB3D5', # EGDP
  'Luzon' = '#D9D9D9', # EGDP
  'Mindanao' = '#C2D6D6', # EGDP
  'Myanmar' = '#B3B3B3', # EGDP
  'Nepal' = '#A6A6A6', # EGDP
  'Palawan' = '#999999', # EGDP
  'Sahul' = '#8C8C8C', # EGDP
  'Sulawesi' = '#808080', # EGDP
  'Vietnam' = '#737373' # EGDP
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
  'Bajau' = 15,  # 2018Cell
  'Saluan' = 15,  # 2018Cell
  'Ami' = 18, # EastSouthAsia
  'Atayal' = 18, # EastSouthAsia
  'Birhor' = 18, # EastSouthAsia
  'Blang' = 18, # EastSouthAsia
  'Blang-DY' = 18, # EastSouthAsia
  'Bougainville' = 18, # EastSouthAsia
  'Cambodian' = 18, # EastSouthAsia
  'Dai' = 18, # EastSouthAsia
  'Kharia' = 18, # EastSouthAsia
  'NAN_Melanesian' = 18, # EastSouthAsia
  'Peranakan_Caucasian' = 18, # EastSouthAsia
  'Peranakan_Chinese' = 18, # EastSouthAsia
  'Peranakan_Chinese_Peranakan_Eurasian' = 18, # EastSouthAsia
  'Peranakan_Eurasian' = 18, # EastSouthAsia
  'Peranakan_Indian' = 18, # EastSouthAsia
  'Thai' = 18, # EastSouthAsia
  'Bangladesh' = 16, # EGDP
  'Borneo' = 16, # EGDP
  'Cebu' = 16, # EGDP
  'India' = 16, # EGDP
  'Luzon' = 16, # EGDP
  'Mindanao' = 16, # EGDP
  'Myanmar' = 16, # EGDP
  'Nepal' = 16, # EGDP
  'Palawan' = 16, # EGDP
  'Sahul' = 16, # EGDP
  'Sulawesi' = 16, # EGDP
  'Vietnam' = 16 # EGDP
)

data_info$Color <- unlist(lapply(data_info$POP, function(sp) {
  regionpop <- data_info$Superpop[data_info$POP == sp][1]
  region_colors[regionpop]
}))

head(data_info)


# Define the marker mapping
# shape_mapping <- setNames(data_info$Shape, data_info$POP)
# # Define the color mapping
# color_mapping <- setNames(data_info$Color, data_info$POP)


# Read PCA Results
pca_data <- read.table("/path/to/malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.snp.geno005.maf005.ld50k.PCA_out.eigenvec", header=FALSE)
# Name columns (Usually PC1, PC2 Column contains V3 And V4)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
eigenvalues <- read.table("/path/to/malay_ho1926_SG10K_PanMLS_BajauSaluan.EGDP.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.snp.geno005.maf005.ld50k.PCA_out.eigenval",header=FALSE)
total_variance <- sum(eigenvalues$V1)
explained_variance <- eigenvalues$V1 / total_variance
explained_variance
# [1] 0.17401735 0.11725917 0.10639215 0.09878987 0.06834546 0.04364081
# [7] 0.04292976 0.03611212 0.03395190 0.03134281 0.02779547 0.02687886
# [13] 0.02666508 0.02616493 0.02492199 0.02344640 0.02336803 0.02289774
# [19] 0.02261303 0.02246707
# Read PCA Data

# Merge PCA data and data_info
merged_data <- merge(pca_data, data_info, by.x = "IID", by.y = "ID", all.x = TRUE)
# Inspect merged results
head(merged_data)


# Plot PCA Plot (PC1 vs PC2)

pdf("malay_ho1926_SG10K_PanMLS_BajauSaluan.chr1_22.merge.Malaysia_sample.90_2018Cell.EGDP.ESA.snp.geno005.maf005.ld50k.PCA_out.geno005.maf005.ld50k.pca.custom.pdf", width=10, height=6)
ggplot(merged_data, aes(x=PC1, y=PC2, color=POP, shape = POP)) +
  geom_point(size=3) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC1 (17.40%)") + ylab("PC2 (11.73%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = region_colors) +  #  Apply custom colors
  scale_shape_manual(values = pop_shapes) +  #  Apply custom markers
  coord_fixed(ratio = 0.5) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )
dev.off()

library(cowplot)
pdf("malay_ho1926_SG10K_PanMLS_BajauSaluan.chr1_22.merge.Malaysia_sample.90_2018Cell.snp.geno005.maf005.ld50k.PCA_out.geno005.maf005.ld50k.pca.pdf", width=6, height=6)
ggplot(merged_data, aes(x = PC1, y = PC2, color = POP, shape = POP)) +
  geom_point(size = 3) +
  xlab("PC1 (18.13%)") + ylab("PC2 (12.24%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = region_colors) +  #  Apply custom colors
  scale_shape_manual(values = pop_shapes) +  #  Apply custom markers
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

