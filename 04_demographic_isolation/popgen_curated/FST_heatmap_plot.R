# Summarize population variation and demographic history.

library(pheatmap)
library(ggplot2)

setwd("/path/to/original_plot")
mat <- read.csv("948.mean.fst", sep = "\t", row.names = 1, check.names = FALSE)

# Calculate distances and cluster
dist_rows <- dist(mat, method = "euclidean")
dist_cols <- dist(t(mat), method = "euclidean")

clust_rows <- hclust(dist_rows, method = "ward.D2")
clust_cols <- hclust(dist_cols, method = "ward.D2")

# Get cluster order
row_order <- clust_rows$order
col_order <- clust_cols$order

# Use Yoruba Move to the beginning (Row)
row_names <- rownames(mat)[row_order]
row_names <- c("Yoruba", setdiff(row_names, "Yoruba"))

# Use Yoruba Move to the beginning (Column)
col_names <- colnames(mat)[col_order]
col_names <- c("Yoruba", setdiff(col_names, "Yoruba"))

# Reorder the matrix
mat_reordered <- mat[row_names, col_names]

color <- colorRampPalette(c("#FEE0D2", "#56B1F7","#132B43"))(100)  #  Blue-to-pink
# Plot
pdf("population_distance_heatmap_rootedYRI.pdf", width = 10, height = 10)
pheatmap(mat_reordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = color,
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Pairwise Population Distance Heatmap",
         cellwidth = 8,     #  Set cell width
         cellheight = 8     #  Set cell height ( Ensure 1:1)
         )
dev.off()
