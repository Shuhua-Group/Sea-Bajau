# Plot statistics in a specified interval, Plot within this interval: Fst, iHS, XPEHH, TajimaD
library(ggplot2)
library(dplyr)

chr<-"chr19"
start<-"10013483"
end<-"10022279"

# corresponding chm13 Coordinates: chr3:96419607-111577885


setwd("/path/to/04.analysis/")


################# 1. Local Fst #################
malay_fst <- read.table("Fst/Bajau_vs_Malay.perSNP.fst.weir.fst", header = TRUE, sep = "\t")
IndoBajau_fst <- read.table("Fst/Bajau_vs_IndoBajau.perSNP.fst.weir.fst", header = TRUE, sep = "\t")

malay_fst$CHR <- gsub("chr", "", malay_fst$CHROM)
malay_fst$CHR <- as.numeric(malay_fst$CHR)
malay_fst_local <- malay_fst %>% filter(CHROM == chr) %>%
  filter(POS >= as.numeric(start)) %>% filter(POS <= as.numeric(end))

IndoBajau_fst$CHR <- gsub("chr", "", IndoBajau_fst$CHROM)
IndoBajau_fst$CHR <- as.numeric(IndoBajau_fst$CHR)
IndoBajau_fst_local <- IndoBajau_fst %>% filter(CHROM == chr) %>%
  filter(POS >= as.numeric(start)) %>% filter(POS <= as.numeric(end))


pdf(paste0("Fst/",chr,"_",start,"_",end,".Fst.Bajau_vs_Malay.perSNP.plot.pdf"),width = 6,height = 3)
gg_manhattan <- ggplot(malay_fst_local, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = factor(CHR %% 2))) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("skyblue3", "darkblue")) +
  geom_vline(xintercept = 10018876, linetype = "dashed", color = "darkred") +
  labs(
    title = "Fst Manhattan Plot",
    x = "Chromosomal Position",
    y = "Fst"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0.1, "lines")
  )
print(gg_manhattan)
dev.off()

pdf(paste0("Fst/",chr,"_",start,"_",end,".Fst.Bajau_vs_IndoBajau.perSNP.plot.pdf"),width = 6,height = 3)
gg_manhattan <- ggplot(IndoBajau_fst_local, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = factor(CHR %% 2))) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("skyblue3", "darkblue")) +
  geom_vline(xintercept = 10018876, linetype = "dashed", color = "darkred") +
  labs(
    title = "Fst Manhattan Plot",
    x = "Chromosomal Position",
    y = "Fst"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0.1, "lines")
  )
print(gg_manhattan)
dev.off()

################# 2. Local iHS #################
ihs_file <- paste0("/path/to/36Bajau_like.GRCh38.dv161.",chr,".ihs.out.norm")
ihs_chr<-read.table(ihs_file, header = FALSE, stringsAsFactors = FALSE)
names(ihs_chr) <- c("CHR", "BP", "ID", "norm_iHS")
ihs_chr$iHS <- as.numeric(ihs_chr$norm_iHS)
ihs_local <- ihs_chr %>%
  filter(BP >= as.numeric(start)) %>% filter(BP <= as.numeric(end))

pdf(paste0("iHS/",chr,"_",start,"_",end,"iHS.Bajau.plot.pdf"),width = 6,height = 3)
# Set colors and threshold lines
gg_manhattan <- ggplot(ihs_local, aes(x = BP, y = norm_iHS)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 10018876, linetype = "dashed", color = "darkred") +
  labs(
    title = "iHS Manhattan Plot",
    x = "Chromosomal Position",
    y = "iHS"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0.1, "lines")
  )
print(gg_manhattan)
dev.off()


################# 3. Local XPEHH #################
malay_XPEHH_data <- read.table(paste0("XPEHH/Malay_XPEHH/36Bajau_vs_37Malay.",chr,".xpehh_output.xpehh.out.norm"), header = TRUE, stringsAsFactors = FALSE)
IndoBajau_XPEHH_data <- read.table(paste0("XPEHH/IndoBajau_XPEHH/36Bajau_vs_IndoBajau.",chr,".xpehh_output.xpehh.out.norm"), header = TRUE, stringsAsFactors = FALSE)
malay_XPEHH_data$pos <- as.numeric(malay_XPEHH_data$pos)
malay_XPEHH_data$normxpehh <- as.numeric(malay_XPEHH_data$normxpehh)
IndoBajau_XPEHH_data$pos <- as.numeric(IndoBajau_XPEHH_data$pos)
IndoBajau_XPEHH_data$normxpehh <- as.numeric(IndoBajau_XPEHH_data$normxpehh)

malay_XPEHH_local <- malay_XPEHH_data %>%
  filter(pos >= as.numeric(start)) %>% filter(pos <= as.numeric(end))
IndoBajau_XPEHH_local <- IndoBajau_XPEHH_data %>%
  filter(pos >= as.numeric(start)) %>% filter(pos <= as.numeric(end))

pdf(paste0("XPEHH/",chr,"_",start,"_",end,".Bajau_vs_IndoBajau.XPEHH.plot.pdf"),width = 6,height = 3)
# Set colors and threshold lines
gg_manhattan <- ggplot(IndoBajau_XPEHH_local, aes(x = pos, y = normxpehh)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 10018876, linetype = "dashed", color = "darkred") +
  labs(
    title = "XPEHH Manhattan Plot",
    x = "Chromosomal Position",
    y = "XP-EHH"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0.1, "lines")
  )
print(gg_manhattan)
dev.off()

pdf(paste0("XPEHH/",chr,"_",start,"_",end,".Bajau_vs_Malay.XPEHH.plot.pdf"),width = 6,height = 3)
# Set colors and threshold lines
gg_manhattan <- ggplot(malay_XPEHH_local, aes(x = pos, y = normxpehh)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 10018876, linetype = "dashed", color = "darkred") +
  labs(
    title = "XPEHH Manhattan Plot",
    x = "Chromosomal Position",
    y = "XP-EHH"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0.1, "lines")
  )
print(gg_manhattan)
dev.off()



################# 4. TajimaD #################
bajau_TajimaD <- read.table("TajimaD/Bajau_50kwindow_TajimaD.Tajima.D", header = TRUE, sep = "\t")

bajau_TajimaD_local <- bajau_TajimaD %>% filter(CHROM == chr) %>%
  # filter(N_SNPS >= 10) %>%
  filter(BIN_START >= as.numeric(start)) %>% filter(BIN_START <= as.numeric(end))

pdf(paste0("TajimaD/",chr,"_",start,"_",end,".TajimaD.Bajau.50k_win.plot.pdf"),width = 6,height = 3)
gg_manhattan <- ggplot(bajau_TajimaD_local, aes(x = BIN_START, y = TajimaD)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = 10018876, linetype = "dashed", color = "darkred") +
  labs(
    title = "TajimaD Manhattan Plot",
    x = "Chromosomal Position",
    y = "TajimaD"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0.1, "lines")
  )
print(gg_manhattan)
dev.off()

################# 5. CMS #################

CMS <- read.csv("/path/to/merged_all_clean.txt",sep="\t",header=TRUE)
CMS_data <- CMS[,c("CHR","POS","ID","CMS")]

CMS_data$CHR2 <- gsub("chr", "", CMS_data$CHR)
CMS_data$CHR2 <- as.numeric(CMS_data$CHR2)

bajau_CMS_local <- CMS_data %>% filter(CHR == chr) %>%
  filter(POS >= as.numeric(start)) %>% filter(POS <= as.numeric(end))

pdf(paste0("CMS/",chr,"_",start,"_",end,".CMS.Bajau.pdf"),width = 6,height = 3)
gg_manhattan <- ggplot(bajau_CMS_local, aes(x = POS, y = CMS)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = 10018876, linetype = "dashed", color = "darkred") +
  labs(
    title = "CMS Manhattan Plot",
    x = "Chromosomal Position",
    y = "CMS"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0.1, "lines")
  )
print(gg_manhattan)
dev.off()

