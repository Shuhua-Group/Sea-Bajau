#conda activate r4
library(ggplot2)
library(RColorBrewer)
library(scales)
library(colorspace)
library(dplyr)

# 1. 102 Malaysia Long-read data [ Keep all records ]
# 2. ho_1926 Short-read data [ Keep selected East and Southeast Asian samples ]
# 3. PanMLS Short-read data [../../reference/PanMLS461/PanMLS461.select.list]
# 4. EGDP Short-read data [ A single level ../../reference/EGDP/sample_info.txt]
# 5. 2018cell Low-coverage short-read data [ Keep all records ]
# 6. Pacific Array data [ Keep all records initially ]
# 7. SG10K Singapore data [ A single level Peranakan_79.Indian_50]


# Read set1:Malay Sample metadata
malay_data_info <- read.table("/path/to/118Sample.info.txt", header=TRUE)
colnames(malay_data_info) <- c("Sample", "FID", "Gender","Superpop")
# Read set2:ho_1926 Sample metadata
ho1926_info <- read.table("/path/to/ho_1926.kp.info.txt", header=FALSE)
colnames(ho1926_info) <- c("Sample", "Subpop", "Superpop")
# Read set3:PanMLS Sample metadata
panmls_info <- read.table("/path/to/PanMLS461.select.sample_info.txt", header=FALSE,sep="\t")
colnames(panmls_info) <- c("Sample", "Superpop","Source")
# Read set4:EGDP Sample metadata
egdp_info <- read.table("/path/to/EGDP.info.txt", header=TRUE)
colnames(egdp_info) <- c("Sample", "Subpop", "Superpop","Source","Language")
# Read set5:2018cell Sample metadata
cell_info <- read.table("/path/to/90_2018Cell.info.txt", header=TRUE)
colnames(cell_info) <- c("Sample", "Subpop", "Superpop","Source","Language")
# Read set6:Pacific Sample metadata
pacific_info <- read.table("/path/to/Pacific.info.txt", header=FALSE)
colnames(pacific_info) <- c("Sample", "Subpop", "Superpop","Source","Language")
# Read set7:SG10K Sample metadata
sg10k_info <- read.table("/path/to/SG10K.Peranakan_79.Indian_50.txt", header=FALSE,sep="\t")
colnames(sg10k_info) <- c("Sample", "Superpop")

# Merge metadata
merged_data_info <- rbind(malay_data_info[,c("Sample", "Superpop")],
                          ho1926_info[,c("Sample", "Superpop")],
                          panmls_info[,c("Sample", "Superpop")],
                          egdp_info[,c("Sample", "Superpop")],
                          cell_info[,c("Sample", "Superpop")],
                          pacific_info[,c("Sample", "Superpop")],
                          sg10k_info[,c("Sample", "Superpop")])

# Get the actual level0 samples
level2_sample <- read.table("/path/to/level2.list", header=FALSE)
colnames(level2_sample) <- c("Sample")
merged_data_info <- merged_data_info[merged_data_info$Sample %in% level2_sample$Sample,]

# Get unique Value
superpops <- unique(merged_data_info$Superpop)

pop_colors <- c(
  'Half-Senoi-Semoq-Beri' = '#FF7F9F', # set1: Malaysia
  'Malay' = '#fc8d59', # set1: Malaysia
  'Negrito-Bateq' = '#d53e4f', # set1: Malaysia
  'Negrito-Mendriq' = '#FF69B4', # set1: Malaysia
  'Quarter-Senoi-Semoq-Beri' = '#FFB6C1', # set1: Malaysia
  'Senoi-Semoq-Beri' = '#b8a3e1', # set1: Malaysia
  'Senoi-Temiar' = '#9e70a9', # set1: Malaysia
  "Senoi-Semai" = '#9370DB', # set1: Malaysia
  'Brunei' = '#fdae61', # set1: Malaysia
  'Bugis' = '#f46d43', # set1: Malaysia
  'Dusun' = '#fca875', # set1: Malaysia
  'Iban' = '#fee08b', # set1: Malaysia
  'Jawa' = '#e58267', # set1: Malaysia
  'Kimaragang' = '#ff9e6a', # set1: Malaysia
  'Sea-Bajau' = '#4169E1', # set1: Malaysia
  "EastAsia" = "#62C15D",    # set2:ho_1926
  "SouthAsia" = "#E78AC3", # set2:ho_1926
  "Africa" = "#D18C7F",     # set2:ho_1926
  "WestEurasia" = "#8DA0CB",   # set2:ho_1926
  "America" = "#B99B60",     # set2:ho_1926
  "Oceania" = "#1aac7c",    # set2:ho_1926
  "CentralAsiaSiberia" = "#1F77B4",      # set2:ho_1926
  "Altaian" = "#62C15D",    # PanMLS
  "Ami" = "#62C15D",        # PanMLS
  "Atayal" = "#62C15D",     # PanMLS
  "Australian" = "#1aac7c", # PanMLS
  "Bateq" = "#fca800",      # PanMLS
  "Birhor" = "#E78AC3",     # PanMLS
  "Bougainville" = "#1aac7c",     # PanMLS
  "Burmese" = "#62C15D",     # PanMLS
  "Cambodian" = "#62C15D",      # PanMLS
  "CheWong" = "#fca800",        # PanMLS
  "Dai" = "#62C15D",     # PanMLS
  "Daur" = "#62C15D", # PanMLS
  "Eskimo_Chaplin" = "#62C15D",      # PanMLS
  "Eskimo_Naukan" = "#62C15D",     # PanMLS
  "Eskimo_Sireniki" = "#62C15D",     # PanMLS
  "Even" = "#62C15D",      # PanMLS
  "Hawaiian" = "#1aac7c",      # PanMLS
  "Hezhen" = "#62C15D",        # PanMLS
  "Igorot" = "#fca800",     # PanMLS
  "Itelman" = "#62C15D", # PanMLS
  "Japanese" = "#62C15D", # PanMLS
  "Jakun" = "#fca800",      # PanMLS
  "Kinh" = "#62C15D",     # PanMLS
  "Korean" = "#62C15D",     # PanMLS
  "Kusunda" = "#fca800",      # PanMLS
  "Kyrgyz" = "#62C15D",        # PanMLS
  "Lahu" = "#62C15D",     # PanMLS
  "Mendriq" = "#fca800", # PanMLS
  "Miao" = "#62C15D",      # PanMLS
  "Mongola" = "#62C15D",     # PanMLS
  "Naxi" = "#62C15D",     # PanMLS
  "Onge" = "#E78AC3",      # PanMLS
  "Papuan" = "#1aac7c",        # PanMLS
  "She" = "#62C15D",     # PanMLS
  "Temiar" = "#fca800", # PanMLS
  "Thai" = "#fca800",      # PanMLS
  "Tibetan" = "#62C15D",     # PanMLS
  "Tu" = "#62C15D",     # PanMLS
  "Tujia" = "#62C15D",      # PanMLS
  "Ulchi" = "#62C15D",        # PanMLS
  "Xibo" = "#62C15D",     # PanMLS
  "Yakut" =  "#62C15D", # PanMLS
  "Yi" = "#62C15D",      # PanMLS
  "Bajau" = "#fca800",     # 2018Cell
  "Saluan" = "#fca800",     # 2018Cell
  "Indonesia" = "#fca800",      # Pacific
  "Borneo" = "#fca800",        # Pacific
  "Philippines" = "#fca800",     # Pacific
  "Peranakan Chinese" = "#fca800", # SG10K
  "Peranakan Indian" = "#fca800",      # SG10K
  "Peranakan Eurasian" = "#fca800",     # SG10K
  "Peranakan Caucasian" = "#fca800",     # SG10K
  "Peranakan Chinese/Peranakan Eurasian" = "#fca800",      # SG10K
  "Indian" = "#E78AC3"        # SG10K
  )






# Read PCA Results
setwd("/path/to/plink")


setwd("/path/to/plink")
pca_data <- read.table("pcs7panel.chrA.level2.geno001.maf005.ld", header=TRUE)
# Name columns (Usually PC1, PC2 Column contains V3 And V4)
colnames(pca_data) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4")
pve <- read.table("pve7panel.chrA.level2.geno001.maf005.ld",header=FALSE)
pve

# V1
# 1 0.024487220
# 2 0.017516330
# 3 0.014621940
# 4 0.011001140
# 5 0.006836504
# 6 0.006458635
# 7 0.006005603
# 8 0.005078410
# 9 0.004749430
# 10 0.004294689

final_data <- merge(pca_data, merged_data_info, by.x = "IID", by.y = "Sample", all.x = TRUE)
malaysia_populations <- c('Half-Senoi-Semoq-Beri', 'Malay', 'Negrito-Bateq',
                          'Negrito-Mendriq', 'Quarter-Senoi-Semoq-Beri',
                          'Senoi-Semoq-Beri', 'Senoi-Semai','Senoi-Temiar', 'Brunei',
                          'Bugis', 'Dusun', 'Iban', 'Jawa', 'Kimaragang',
                          'Sea-Bajau')
final_data$Group <- ifelse(final_data$Superpop %in% malaysia_populations, "Malaysia", "Other")
final_data <- final_data[!is.na(final_data$Superpop), ]
final_data <- final_data[order(final_data$Group == "Malaysia"), ]



library(cowplot)

pdf("7panel.chrA.level2.flashpca.pdf", width=8, height=8)
ggplot(final_data, aes(x=PC1, y=PC2, fill=Superpop)) +
  geom_point(shape = 21, size = 3, aes(color = Group)) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC1(2.45%)") + ylab("PC2(1.75%)") +
  ggtitle("PCA Plot") +
  scale_fill_manual(values = pop_colors) +  #  Apply custom colors
  scale_color_manual(values = c("Malaysia" = "black", "Other" = NA)) +
  coord_fixed(ratio = 0.715) +  #  Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )+
  theme(legend.position = "none")  #  Hide the legend
dev.off()

pca_plot <- ggplot(final_data, aes(x=PC1, y=PC2, fill=Superpop)) +
  geom_point(shape = 21, size = 3, aes(color = Group)) +
  #geom_text(aes(label=IID), vjust=-0.5, hjust=0.5, size=3) + # Add IID Labels
  xlab("PC1(2.45%)") + ylab("PC2(1.75%)") +
  ggtitle("PCA Plot") +
  scale_fill_manual(values = pop_colors) +  #  Apply custom colors
  scale_color_manual(values = c("Malaysia" = "black", "Other" = NA)) +
  #coord_fixed(ratio = 1) + # Use a fixed aspect ratio
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=16, face="bold"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=10)
  )
legend_plot <- get_legend(pca_plot)

# Save the legend
ggsave("7panel.chrA.level2.flashpca.legend.pdf", legend_plot, width = 10, height = 6)



