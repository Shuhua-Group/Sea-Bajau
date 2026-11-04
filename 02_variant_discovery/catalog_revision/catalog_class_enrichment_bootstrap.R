#!/usr/bin/env Rscript
# Block-bootstrap enrichment of each catalog's strictly unique alleles by variant class.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2 || length(args) > 3) stop("Usage: catalog_class_enrichment_bootstrap.R detailed.tsv out_prefix [n_boot=1000]")
d <- read.delim(args[1], check.names=FALSE)
out <- args[2]; B <- ifelse(length(args)==3, as.integer(args[3]), 1000L)
catalogs <- c("APC", "CPC2", "HPRC2")
blocks <- split(seq_len(nrow(d)), d$BLOCK_1MB)
classes <- sort(unique(d$CLASS))
calc <- function(ix, catalog, klass) {
  u <- d[[catalog]][ix] == 1 & rowSums(d[ix, catalogs]) == 1
  k <- d$CLASS[ix] == klass
  a <- sum(u & k); b <- sum(u & !k); c <- sum(!u & k); e <- sum(!u & !k)
  log(((a+.5)*(e+.5))/((b+.5)*(c+.5)))
}
set.seed(20260907)
ans <- list(); z <- 1L
for (catalog in catalogs) for (klass in classes) {
  obs <- calc(seq_len(nrow(d)), catalog, klass)
  boot <- replicate(B, { ix <- unlist(sample(blocks, length(blocks), replace=TRUE), use.names=FALSE); calc(ix, catalog, klass) })
  ans[[z]] <- data.frame(CATALOG=catalog, CLASS=klass, logOR=obs, OR=exp(obs),
    OR_low=exp(quantile(boot,.025)), OR_high=exp(quantile(boot,.975))); z <- z+1L
}
r <- do.call(rbind, ans); write.table(r, paste0(out, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
p <- ggplot(r, aes(CLASS, OR, colour=CATALOG)) + geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(aes(ymin=OR_low,ymax=OR_high), position=position_dodge(.55), width=.2) +
  geom_point(position=position_dodge(.55), size=2.4) + scale_y_log10() + theme_classic(base_size=12) +
  labs(x=NULL,y="Odds ratio: strictly catalog-specific vs all other alleles") + theme(axis.text.x=element_text(angle=35,hjust=1))
ggsave(paste0(out, ".pdf"),p,width=8,height=4.8)
