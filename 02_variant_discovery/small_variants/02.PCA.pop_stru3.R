#conda activate r4
library(ggplot2)
library(RColorBrewer)
library(scales)
library(colorspace)
# This stru3 Add PanMLS Population, This version does not distinguish subpopulations, Consider only Superpop
# However PanMLS sample metadata do not contain Superpop, Use the available fine-scale labels


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
  'Half-Senoi-Semoq-Beri' = 19,  #  Triangle
  'Malay' = 19,
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
panmls_info <- read.table("/path/to/sample_info.txt", header=TRUE)
colnames(panmls_info) <- c("Sample", "Subpop", "Source")
# Get unique Value
subpops <- unique(panmls_info$Subpop)
# Assign a base color to each population
set.seed(123) #  Make assignments reproducible
mls_colors <- scales::hue_pal()(length(subpops)) #  Use different hues

mls_shapes <- (0:(length(subpops) - 1)) %% 19 #  Generate  0  To  18  repeating marker sequence
subpop_mapping <- data.frame(
  Subpop = subpops,
  Color = mls_colors,
  Shape = mls_shapes
)

# Convert the format to IID FID Gender Population
# Create a data frame
panmls_transformed <- panmls_info[, c("Sample", "Sample", "Subpop")]
panmls_transformed$Gender <- NA  #  Set  Gender  columns:  NA
# Rename columns
colnames(panmls_transformed) <- c("IID", "FID", "Population", "Gender")
panmls_transformed <- panmls_transformed[, colnames(malay_data_info)]
head(panmls_transformed)

# Read PCA Results
pca_data <- read.table("/path/to/malay_panmls.snp.geno005.maf005.ld50k.PCA_out.eigenvec", header=FALSE)
# Name columns (Usually PC1, PC2 Column contains V3 And V4)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4")
eigenvalues <- read.table("/path/to/malay_panmls.snp.geno005.maf005.ld50k.PCA_out.eigenval",header=FALSE)
total_variance <- sum(eigenvalues$V1)
explained_variance <- eigenvalues$V1 / total_variance

# [1] 0.24947346 0.14077465 0.11031688 0.08343741 0.07457347 0.05559221
# [7] 0.04013167 0.03100730 0.02908085 0.02133092 0.02073966 0.02039019
# [13] 0.01852637 0.01803486 0.01580952 0.01540776 0.01437479 0.01427498
# [19] 0.01362732 0.01309571

# Merge sample metadata with PCA Results
combined_data <- rbind(panmls_transformed,malay_data_info)
head(combined_data)
combined_data["FID"] = 0
merged_data <- merge(pca_data, combined_data, by=c("FID", "IID"))


# Join plot colors and markers
malay_colors_shapes <- data.frame(Subpop = names(malay_colors), Color = malay_colors,Shape = malay_shapes)
combined_colors_shapes <- rbind(subpop_mapping, malay_colors_shapes)
color_mapping <- setNames(combined_colors_shapes$Color, combined_colors_shapes$Subpop)
shape_mapping <- setNames(combined_colors_shapes$Shape, combined_colors_shapes$Subpop)

# Plot PCA Plot (PC1 vs PC2)

# Hide the legend
pdf("/path/to/118Malay.panmls.geno005.maf005.ld50k.pca.superpop.pdf", width=8, height=8)
ggplot(merged_data, aes(x=PC1, y=PC2, color=Population, shape = Population)) +
  geom_point(size=3) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC1(24.95%)") + ylab("PC2(14.08%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = color_mapping) +  #  Apply custom colors
  scale_shape_manual(values = shape_mapping) +  #  Apply custom markers
  coord_fixed(ratio = 1) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )+
  theme(legend.position = "none")  #  Hide the legend
dev.off()

# Draw the legend separately
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
ggsave("118Malay.panmls.legend_plot.pdf", legend_plot, width = 15, height = 6)
