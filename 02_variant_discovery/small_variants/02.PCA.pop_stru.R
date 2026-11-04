#conda activate r4
library(ggplot2)
library(RColorBrewer)
library(scales)
library(colorspace)

# Read 118Malay Sample metadata
malay_data_info <- read.table("/path/to/118Sample.info.txt", header=TRUE)
malay_colors <- c(
  'Half-Senoi-Semoq-Beri' = '#FF7F9F',
  'Malay' = '#4169E1',
  'Negrito-Bateq' = '#d53e4f',
  'Negrito-Mendriq' = '#FF69B4',
  'Quarter-Senoi-Semoq-Beri' = '#FFB6C1',
  'Senoi-Semai' = '#9370DB',
  'Senoi-Semoq-Beri' = '#b8a3e1',
  'Senoi-Temiar' = '#9e70a9',
  'Brunei' = '#007BFF',
  'Bugis' = '#1E90FF',
  'Dusun' = '#6495ED',
  'Iban' = '#ADD8E6',
  'Jawa' = '#00BFFF',
  'Kimaragang' = '#87CEEB',
  'Sea-Bajau' = '#fc8d59'
)

malay_shapes <- c(
  'Half-Senoi-Semoq-Beri' = 17,  #  Triangle
  'Malay' = 19,  #  Filled circle
  'Negrito-Bateq' = 15,  #  Filled square
  'Negrito-Mendriq' = 18,  #  Filled diamond
  'Quarter-Senoi-Semoq-Beri' = 1,  #  Open circle
  'Senoi-Semai' = 2,  #  Open triangle
  'Senoi-Semoq-Beri' = 0,  #  Open square
  'Senoi-Temiar' = 5,  #  Open pentagon
  'Brunei' = 8,  #  Star
  'Bugis' = 3,  #  Plus sign
  'Dusun' = 4,  #  Cross
  'Iban' = 6,  #  Hexagon
  'Jawa' = 7,  #  Six-pointed star
  'Kimaragang' = 12,  #  Thick cross
  'Sea-Bajau' = 10  #  Filled dot
)

# Read ho_1926 Sample metadata
ho1926_info <- read.table("/path/to/ho_1926.kp.info.txt", header=FALSE)
colnames(ho1926_info) <- c("Sample", "Subpop", "Superpop")
# Get unique Value
superpops <- unique(ho1926_info$Superpop)
subpops <- unique(ho1926_info$Subpop)

# Assign a base color to each population
# Set colors explicitly
superpop_colors <- c(
  "EastAsia" = "#FFCC60",    #  Orange
  "CentralAsiaSiberia" = "#1F77B4",      #  Blue
  "Africa" = "#D18C7F",     #  Green
  "SouthAsia" = "#99C2B7", #  Red
  "Oceania" = "#62C15D",    #  Purple
  "America" = "#B99B60",     #  Brown
  "WestEurasia" = "#F3A9D3"     #  Pink-purple
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
pca_data <- read.table("/path/to/malay_ho_1926.snp.geno005.maf005.ld50k.PCA_out.eigenvec", header=FALSE)
# Name columns (Usually PC1, PC2 Column contains V3 And V4)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4")
eigenvalues <- read.table("/path/to/malay_ho_1926.snp.geno005.maf005.ld50k.PCA_out.eigenval",header=FALSE)
total_variance <- sum(eigenvalues$V1)
explained_variance <- eigenvalues$V1 / total_variance

# [1] 0.450681719 0.276522658 0.045916797 0.040256032 0.031134207 0.022233592
# [7] 0.018940724 0.017864852 0.017339233 0.012582592 0.011321611 0.007965708
# [13] 0.007806472 0.006763847 0.006255194 0.006217853 0.005712348 0.005031963
# [19] 0.004805533 0.004647064

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

pdf("118Malay.ho_1926.geno005.maf005.ld50k.pca.custom.pdf", width=20, height=8)
ggplot(merged_data, aes(x=PC1, y=PC2, color=Population, shape = Population)) +
  geom_point(size=3) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC1 (45.07%)") + ylab("PC2 (27.65%)") +
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
pdf("118Malay.ho_1926.geno005.maf005.ld50k.pca.pdf", width=8, height=8)
ggplot(merged_data, aes(x = PC1, y = PC2, color = Population, shape = Population)) +
  geom_point(size = 3) +
  xlab("PC1:45.07%") + ylab("PC2:27.65%") +
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
pca_plot <- ggplot(merged_data, aes(x = PC1, y = PC2, color = Population, shape = Population)) +
  geom_point(size = 3) +
  xlab("PC1:45.07%") + ylab("PC2:27.65%") +
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

ggsave("legend_plot.pdf", legend_plot, width = 20, height = 6)

