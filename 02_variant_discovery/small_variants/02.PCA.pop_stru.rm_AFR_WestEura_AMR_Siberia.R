#conda activate r4
library(ggplot2)
library(RColorBrewer)
library(scales)
library(colorspace)

# Read 118Malay Sample metadata
malay_data_info <- read.table("/path/to/118Sample.info.txt", header=TRUE)
malay_colors <- c(
  'Half-Senoi-Semoq-Beri' = '#FF7F9F',
  'Malay' = '#fc8d59',
  'Negrito-Bateq' = '#d53e4f',
  'Negrito-Mendriq' = '#FF69B4',
  'Quarter-Senoi-Semoq-Beri' = '#FFB6C1',
  'Senoi-Semai' = '#9370DB',
  'Senoi-Semoq-Beri' = '#b8a3e1',
  'Senoi-Temiar' = '#9e70a9',
  'Brunei' = '#fdae61',
  'Bugis' = '#f46d43',
  'Dusun' = '#fca875',
  'Iban' = '#fee08b',
  'Jawa' = '#e58267',
  'Kimaragang' = '#ff9e6a',
  'Sea-Bajau' = '#4169E1'
)


malay_shapes <- c(
  'Half-Senoi-Semoq-Beri' = 19,  #  Triangle
  'Malay' = 19,  #  Filled circle
  'Negrito-Bateq' = 19,  #  Filled square
  'Negrito-Mendriq' = 19,  #  Filled diamond
  'Quarter-Senoi-Semoq-Beri' = 19,  #  Open circle
  'Senoi-Semai' = 19,  #  Open triangle
  'Senoi-Semoq-Beri' = 19,  #  Open square
  'Senoi-Temiar' = 19,  #  Open pentagon
  'Brunei' = 19,  #  Star
  'Bugis' = 19,  #  Plus sign
  'Dusun' = 19,  #  Cross
  'Iban' = 19,  #  Hexagon
  'Jawa' = 19,  #  Six-pointed star
  'Kimaragang' = 19,  #  Thick cross
  'Sea-Bajau' = 19  #  Filled dot
)

# Read ho_1926 Sample metadata
ho1926_info <- read.table("/path/to/ho_1926.kp.info.rm_AFR_WestEura_AMR_Siberia.txt", header=FALSE)
colnames(ho1926_info) <- c("Sample", "Subpop", "Superpop")
# Get unique Value
superpops <- unique(ho1926_info$Superpop)
subpops <- unique(ho1926_info$Subpop)

# Assign a base color to each population
# Set colors explicitly
superpop_colors <- c(
  "EastAsia" = "#62C15D",    #  Orange
  "SouthAsia" = "#E78AC3", #  Pale green
  "Oceania" = "#FFCC60"    #  Green
)
# ho1926_colors <- brewer.pal(length(ho1926_suppop), "Set2")
# "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
names(superpop_colors) <- superpops

population_colors <- list()
for (superpop in superpops) {
  # Get the current Superpop subpopulations
  subpop_in_group <- unique(ho1926_info$Subpop[ho1926_info$Superpop == superpop])
  # Get the base color
  base_color <- superpop_colors[superpop]
  # Generate distinct colors (Light to dark)
  subpop_colors <- colorRampPalette(c(lighten(base_color, 0.4), base_color, darken(base_color, 0.4)))(length(subpop_in_group))
  # For the current Superpop assign colors to subpopulations
  names(subpop_colors) <- subpop_in_group
  population_colors[[superpop]] <- subpop_colors
}
# Map the generated colors to the data
ho1926_info$Color <- unlist(lapply(ho1926_info$Subpop, function(sp) {
  superpop <- ho1926_info$Superpop[ho1926_info$Subpop == sp][1]
  population_colors[[superpop]][sp]
}))
head(ho1926_info)


# Assign a base marker to each population
shapes <- c(16, 17, 15, 18, 19, 20, 8, 3, 7, 4, 9, 10, 11, 12, 13, 14, 5, 6)
set.seed(42)  #  Set the random seed for reproducibility
subpop_shapes <- sample(shapes, length(subpops), replace = TRUE)
subpop_shape_map <- setNames(subpop_shapes, unique(ho1926_info$Subpop))
# Map markers to ho1926_info
ho1926_info$Shape <- subpop_shape_map[ho1926_info$Subpop]

# Convert the format to IID FID Gender Population
# Create a data frame
ho1926_transformed <- ho1926_info[, c("Sample", "Sample", "Subpop")]
ho1926_transformed$Gender <- NA  #  Set  Gender  columns:  NA
# Rename columns
colnames(ho1926_transformed) <- c("IID", "FID", "Population", "Gender")
ho1926_transformed <- ho1926_transformed[, colnames(malay_data_info)]
head(ho1926_transformed)

# Read PCA Results
pca_data <- read.table("/path/to/malay_ho_1926.rm_AFR_WestEura_AMR_Siberia.snp.geno005.maf005.ld50k.PCA_out.eigenvec", header=FALSE)
# Name columns (Usually PC1, PC2 Column contains V3 And V4)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4")
eigenvalues <- read.table("/path/to/malay_ho_1926.rm_AFR_WestEura_AMR_Siberia.snp.geno005.maf005.ld50k.PCA_out.eigenval",header=FALSE)
total_variance <- sum(eigenvalues$V1)
explained_variance <- eigenvalues$V1 / total_variance
explained_variance
# [1] 0.43423474 0.10925537 0.07822567 0.07248692 0.04750945 0.03402738
# [7] 0.02460880 0.02158128 0.01892283 0.01837068 0.01756245 0.01668111
# [13] 0.01638790 0.01538824 0.01498956 0.01275900 0.01205700 0.01190341
# [19] 0.01167764 0.01137059

# Merge sample metadata with PCA Results
combined_data <- rbind(ho1926_transformed,malay_data_info)
head(combined_data)
combined_data["FID"] = 0
merged_data <- merge(pca_data, combined_data, by=c("FID", "IID"))

# Join plot colors and markers
ho1926_colors_shapes <- ho1926_info[, c("Subpop", "Color", "Shape")]
malay_colors_shapes <- data.frame(Subpop = names(malay_colors), Color = malay_colors,Shape = malay_shapes)
combined_colors_shapes <- rbind(ho1926_colors_shapes, malay_colors_shapes)
color_mapping <- setNames(combined_colors_shapes$Color, combined_colors_shapes$Subpop)
shape_mapping <- setNames(combined_colors_shapes$Shape, combined_colors_shapes$Subpop)


# Plot PCA Plot (PC1 vs PC2)

pdf("118Malay.ho_1926_rm_AFR_WestEura_AMR_Siberia.geno005.maf005.ld50k.pca.custom.pdf", width=20, height=8)
ggplot(merged_data, aes(x=PC1, y=PC2, color=Population, shape = Population)) +
  geom_point(size=3) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC1 (43.42%)") + ylab("PC2 (10.93%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = color_mapping) +  #  Apply custom colors
  scale_shape_manual(values = shape_mapping) +  #  Apply custom markers
  coord_fixed(ratio = 0.5) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )
dev.off()

library(cowplot)
pdf("118Malay.ho_1926_rm_AFR_WestEura_AMR_Siberia.geno005.maf005.ld50k.pca.pdf", width=8, height=8)
ggplot(merged_data, aes(x = PC1, y = PC2, color = Population, shape = Population)) +
  geom_point(size = 3) +
  xlab("PC1 (43.42%)") + ylab("PC2 (10.93%)") +
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
pca_plot <- ggplot(merged_data, aes(x = PC1, y = PC2, color = Population, shape = Population)) +
  geom_point(size = 3) +
  xlab("PC1 (43.42%)") + ylab("PC2 (10.93%)") +
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

ggsave("118Malay.ho_1926_rm_AFR_WestEura_AMR_Siberia.geno005.maf005.ld50k.pca.legend_plot.pdf", legend_plot, width = 10, height = 6)

