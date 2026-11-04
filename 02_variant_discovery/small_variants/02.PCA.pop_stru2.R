#conda activate r4
library(ggplot2)
library(RColorBrewer)
library(scales)
library(colorspace)
# This stru2 This version ignores ho_1926 subpopulations, Consider only Superpop

# Read 118Malay Sample metadata
malay_data_info <- read.table("/path/to/118Sample.info.txt", header=TRUE)
pop_colors <- c(
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
  'Sea-Bajau' = '#4169E1',
  "EastAsia" = "#62C15D",    #  Orange
  "CentralAsiaSiberia" = "#1F77B4",      #  Blue
  "Africa" = "#D18C7F",     #  Green
  "SouthAsia" = "#E78AC3", #  Red
  "Oceania" = "#FFCC60",    #  Purple
  "America" = "#B99B60",     #  Brown
  "WestEurasia" = "#8DA0CB"    #  Pink-purple
)



pop_shapes <- c(
  'Half-Senoi-Semoq-Beri' = 17,  #  Triangle
  'Malay' = 9,
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
  'Sea-Bajau' = 10,  #  Filled dot
  "EastAsia" = 20,
  "CentralAsiaSiberia" = 20,
  "Africa" = 20,
  "SouthAsia" = 20,
  "Oceania" = 20,
  "America" = 20,
  "WestEurasia" = 20
)

# Read ho_1926 Sample metadata
ho1926_info <- read.table("/path/to/ho_1926.kp.info.txt", header=FALSE)
colnames(ho1926_info) <- c("Sample", "Subpop", "Superpop")
# Get unique Value
superpops <- unique(ho1926_info$Superpop)
subpops <- unique(ho1926_info$Subpop)

# Assign a base color to each population


# Convert the format to IID FID Gender Population
# Create a data frame
ho1926_transformed <- ho1926_info[, c("Sample", "Sample", "Superpop")]
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
explained_variance

# [1] 0.450681719 0.276522658 0.045916797 0.040256032 0.031134207 0.022233592
# [7] 0.018940724 0.017864852 0.017339233 0.012582592 0.011321611 0.007965708
# [13] 0.007806472 0.006763847 0.006255194 0.006217853 0.005712348 0.005031963
# [19] 0.004805533 0.004647064

# Merge sample metadata with PCA Results
combined_data <- rbind(ho1926_transformed,malay_data_info)
head(combined_data)
combined_data["FID"] = 0
merged_data <- merge(pca_data, combined_data, by=c("FID", "IID"))

# Plot PCA Plot (PC1 vs PC2)

pdf("/path/to/118Malay.ho_1926.geno005.maf005.ld50k.pca.superpop.pdf", width=12, height=8)
ggplot(merged_data, aes(x=PC1, y=PC2, color=Population, shape = Population)) +
  geom_point(size=3) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC1(45.07%)") + ylab("PC2(27.65%)") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = pop_colors) +  #  Apply custom colors
  scale_shape_manual(values = pop_shapes) +  #  Apply custom markers
  coord_fixed(ratio = 1) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )
dev.off()

