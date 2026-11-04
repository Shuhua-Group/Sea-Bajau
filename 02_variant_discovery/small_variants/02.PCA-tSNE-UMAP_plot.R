#conda activate r4
library(ggplot2)

population_colors <- c(
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
#FF7F9F
#4169E1
#fc8d59
#66C2A5
#D4B9DA
#8DA0CB
#E78AC3
#A6D854


population_shapes <- c(
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

# population_colors <- c(
# 'Half-Senoi-Semoq-Beri' = '#5e3c99', # Dark purple
# 'Malay' = '#fc8d59', # Orange
# 'Negrito-Bateq' = '#d53e4f', # Keep red unchanged
# 'Negrito-Mendriq' = '#FF69B4' # Keep light pink unchanged
# 'Quarter-Senoi-Semoq-Beri' = '#9370DB', # Purple (Bright)
# 'Senoi-Semai' = '#7B68EE', # Violet blue
# 'Senoi-Semoq-Beri' = '#BA55D3', # Lilac
# 'Senoi-Temiar' = '#DA70D6', # Light magenta
# 'Brunei' = '#FFD700', # Golden yellow
# 'Bugis' = '#FFA07A', # Light orange
# 'Dusun' = '#FF8C00', # Dark orange
# 'Iban' = '#FF6347', # Tomato red
# 'Jawa' = '#FFDAB9', # Light peach
# 'Kimaragang' = '#F4A460', # Sandy brown
# 'Sea-Bajau' = '#4169E1' # Blue
#)


# Read PCA Results
#pca_data <- read.table("/path/to/67West_51East.chm13v2.dv161.chr1_22.snp.geno005.maf005.ld50k.PCA_out.eigenvec", header=FALSE)
pca_data <- read.table("/path/to/67West_51East.rm16.geno005.maf005.ld50k.PCA_out.eigenvec", header=FALSE)

population_data <- read.table("/path/to/118Sample.info.txt", header=TRUE)
# Name columns (Usually PC1, PC2 Column contains V3 And V4)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4")
eigenvalues <- read.table("/path/to/67West_51East.rm16.geno005.maf005.ld50k.PCA_out.eigenval",header=FALSE)
total_variance <- sum(eigenvalues$V1)
explained_variance <- eigenvalues$V1 / total_variance
explained_variance

### Calculate the variance explained by each principal component
## Before removal 16 first-degree relatives
# [1] 0.27649636 0.10552549 0.06116086 0.05942854 0.04244740 0.04052418
# [7] 0.03862323 0.03403197 0.03363327 0.03301068 0.03170119 0.03117879
# [13] 0.02956440 0.02783394 0.02700973 0.02693263 0.02624366 0.02563015
# [19] 0.02498880 0.02403472

## After removal 16 first-degree relatives
# [1] 0.26919481 0.11636421 0.04856485 0.04608390 0.04199713 0.03861069
# [7] 0.03720036 0.03510753 0.03425124 0.03296553 0.03185489 0.03147848
# [13] 0.03112602 0.03074146 0.02995448 0.02988924 0.02931099 0.02876204
# [19] 0.02834434 0.02819780

merged_data <- merge(pca_data, population_data, by=c("FID", "IID"))

# Plot PCA Plot (PC1 vs PC2)
#tSNE Plot
library(Rtsne)
library(umap)
library(ggplot2)
data <- read.table("102Samples.geno005.maf005.ld50k.GT.raw",header=TRUE)
genotype_data <- as.matrix(data[,-c(1,2,3,4,5,6)])
set.seed(42)

# Impute missing values with column means _result$Y)
for (i in 1:ncol(genotype_data)) {
  genotype_data[is.na(genotype_data[, i]), i] <- mean(genotype_data[, i], na.rm = TRUE)
}
tsne_result <- Rtsne(genotype_data, perplexity=2)
tsne_df <- as.data.frame(tsne_result$Y)
colnames(tsne_df) <- c("Dim1", "Dim2")
sample_names <- data[, 1]
tsne_df$Sample <- sample_names


population_data <- read.table("/path/to/118Sample.info.txt", header=TRUE)
tsne_merged <- merge(tsne_df, population_data, by.x = "Sample", by.y = "IID", all.x = TRUE)
pdf("118Samples.geno005.maf005.ld50k.TSNE.custom.pdf", width=8, height=8)
ggplot(tsne_merged, aes(x=Dim1, y=Dim2, color=Population)) +
  geom_point(size=3) +
  geom_text(aes(label=Sample), vjust=-0.5, hjust=0.5, size=3) +  #  Add  Sample  Labels
  xlab("t-SNE 1") + ylab("t-SNE 2") +
  ggtitle("t-SNE Plot") +
  scale_color_manual(values = population_colors) +  #  Apply custom colors
  coord_fixed(ratio = 1) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )
dev.off()

##UMAP Plot

umap_result <- umap(genotype_data)
umap_df <- data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- data$IID
umap_merged <- merge(umap_df, population_data,by.x = "Sample", by.y = "IID", all.x = TRUE)
pdf("118Samples.geno005.maf005.ld50k.UMAP.custom.pdf", width=8, height=8)
ggplot(umap_merged, aes(x=UMAP1, y=UMAP2, color=Population)) +
  geom_point(size=3) +
  geom_text(aes(label=Sample), vjust=-0.5, hjust=0.5, size=3) +  #  Add  IID  Labels
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("UMAP Plot") +
  scale_color_manual(values = population_colors) +  #  Apply custom colors
  coord_fixed(ratio = 1) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )
dev.off()







## GCTA## Identical results
# /path/to/gcta64 --bfile chr1-22.snp.rm9.mind005.geno005.maf005 --make-grm --out chr1-22.snp.rm9.mind005.geno005.maf005.gcta_grm
# /path/to/gcta64 --grm chr1-22.snp.rm9.mind005.geno005.maf005.gcta_grm --pca 10 --out chr1-22.snp.rm9.mind005.geno005.maf005.gcta_grm.pca

# library(ggplot2)

# # Read PCA Results
# pca_data <- read.table("chr1-22.snp.rm9.mind005.geno005.maf005.gcta_grm.pca.eigenvec", header=FALSE)
# population_data <- read.table("/path/to/76sample.PASS_67.info.txt2", header=TRUE)
# # Name columns (Usually PC1, PC2 Column contains V3 And V4)
# colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4")
# eigenvalues <- read.table("chr1-22.snp.rm9.mind005.geno005.maf005.gcta_grm.pca.eigenval",header=FALSE)
# total_variance <- sum(eigenvalues$V1)
# explained_variance <- eigenvalues$V1 / total_variance

# merged_data <- merge(pca_data, population_data, by=c("FID", "IID"))

# # Plot PCA Plot (PC1 vs PC2)
# pdf("58PASS_merged.mind005.geno005.maf005.ld50k.gcta_pca.pdf",width=12,height=8)
# ggplot(merged_data, aes(x=PC1, y=PC2, color=Population)) +
# geom_point() +
# xlab("PC1") + ylab("PC2") +
# ggtitle("PCA Plot")

# dev.off()
