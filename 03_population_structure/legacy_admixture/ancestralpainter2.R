# Prepare population datasets and summarize genetic structure.

library(AncestryPainterV2)
exp_q <- read.table("/path/to/admixture.5.Q", header = F)
exp_ind <- read.table("/path/to/admixture.info.sorted.txt", stringsAsFactors = F, header = F)
exp_q[201:205, 1:5]
head(exp_ind)
exp_cols <- readLines("/path/to/admixture.colors.txt")
exp_cols_df <- data.frame(colors = exp_cols, stringsAsFactors = FALSE)$colors
print(exp_cols_df)


pdf("/path/to/admixture.K5.pdf", width = 45, height = 45)
exp_tars <- c("Bajau", "Malay", "Negrito")
sectorplot(Q = exp_q, ind = exp_ind, ancescols = exp_cols_df, target = exp_tars)
dev.off()

