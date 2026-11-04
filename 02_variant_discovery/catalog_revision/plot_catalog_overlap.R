#!/usr/bin/env Rscript
# Exact normalized allele overlap: main plot. SV tolerant matching is a separate sensitivity analysis.
suppressPackageStartupMessages({ library(ggplot2); library(VennDiagram); library(grid) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: plot_catalog_overlap.R membership.detailed.tsv out_prefix")
tab <- read.delim(args[1], check.names = FALSE)
out <- args[2]
set_sizes <- c(APC = sum(tab$APC == 1), CPC2 = sum(tab$CPC2 == 1), HPRC2 = sum(tab$HPRC2 == 1))
pairwise <- c(sum(tab$APC == 1 & tab$CPC2 == 1), sum(tab$APC == 1 & tab$HPRC2 == 1), sum(tab$CPC2 == 1 & tab$HPRC2 == 1))
triple <- sum(tab$APC == 1 & tab$CPC2 == 1 & tab$HPRC2 == 1)
pdf(paste0(out, ".venn.pdf"), width = 7, height = 7)
grid.draw(draw.triple.venn(area1=set_sizes["APC"], area2=set_sizes["CPC2"], area3=set_sizes["HPRC2"],
  n12=pairwise[1], n13=pairwise[2], n23=pairwise[3], n123=triple,
  category=c("APC-P1", "CPC2", "HPRC2"), fill=c("#D95F02", "#1B9E77", "#7570B3"), alpha=0.42,
  cex=1.2, cat.cex=1.3, cat.col=c("#D95F02", "#1B9E77", "#7570B3")))
dev.off()

class_tab <- aggregate(rep(1, nrow(tab)), by=list(CLASS=tab$CLASS, MEMBERSHIP=tab$MEMBERSHIP), FUN=sum)
names(class_tab)[3] <- "N"
class_tab$FRACTION_WITHIN_CLASS <- class_tab$N / ave(class_tab$N, class_tab$CLASS, FUN=sum)
p <- ggplot(class_tab, aes(CLASS, FRACTION_WITHIN_CLASS, fill=MEMBERSHIP)) +
  geom_col(position="fill", width=.8) + scale_y_continuous(labels=scales::percent_format()) +
  labs(x=NULL, y="Fraction of normalized alleles", fill="Catalog membership") + theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=35, hjust=1))
ggsave(paste0(out, ".class_membership.pdf"), p, width=8, height=4.8)
