############################################################
# Assembly N50 vs Variant Discovery
#
# Purpose:
# Evaluate whether assembly contiguity (N50) is associated
# with downstream variant discovery.
#
# Reviewer concern:
# Could assembly-quality heterogeneity systematically affect
# variant discovery rates?
#
# Data structure:
#
# Haplotype N50 AssemblySize Sample SNP Indel Gender
# APC0001.hap1 92400346 ... APC0001 5096645 1147119 Male
# APC0001.hap2 82664270 ... APC0001 5096645 1147119 Male
#
# Important:
#
# SNP and Indel:
# sample-level quantities
# duplicated on the two haplotype rows
#
# INS / DEL:
# haplotype-level quantities
#
# INStelo / DELtelo:
# haplotype-level quantities
#
# INScenSat / DELcenSat:
# haplotype-level quantities
#
############################################################


############################################################
# 1. Load packages
############################################################

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)


############################################################
# 2. Parameters
############################################################

input_file <- "/path/to/130APC_n50_genomesize.txt"

setwd("/path/to/asm_qv")
# Output directory

output_prefix <- "APC_N50_variant"

# output directory
outdir <- "N50_analysis"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 3. Read data
############################################################

dat <- fread(
  input_file,
  sep = "\t",
  header = TRUE
)


############################################################
# 4. Basic checks
############################################################

required_columns <- c(
  "Haplotype",
  "N50",
  "AssemblySize",
  "Sample",
  "SNP",
  "Indel",
  "Gender",
  "INS",
  "DEL",
  "DELtelo",
  "INStelo",
  "INScenSat",
  "DELcenSat"
)


missing_columns <- setdiff(
  required_columns,
  colnames(dat)
)


if (length(missing_columns) > 0) {

  stop(
    paste(
      "Missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  )

}


############################################################
# 5. Convert numeric columns
############################################################

numeric_columns <- c(
  "N50",
  "AssemblySize",
  "SNP",
  "Indel",
  "INS",
  "DEL",
  "DELtelo",
  "INStelo",
  "INScenSat",
  "DELcenSat"
)


dat[
  ,
  (numeric_columns) :=
    lapply(
      .SD,
      as.numeric
    ),
  .SDcols = numeric_columns
]


############################################################
# 6. Check number of haplotypes per sample
############################################################

hap_count <- dat[
  ,
  .(
    n_hap = uniqueN(Haplotype)
  ),
  by = Sample
]


if (any(hap_count$n_hap != 2)) {

  warning(
    "Some samples do not have exactly two haplotypes."
  )

}


############################################################
# 7. Haplotype-level analysis
############################################################
#
# These variables are genuinely haplotype-specific:
#
# INS
# DEL
# INStelo
# DELtelo
# INScenSat
# DELcenSat
#
# Therefore:
#
# N50(haplotype)
# vs
# variant count(haplotype)
#
# is statistically appropriate.
#
############################################################


hap_df <- dat


############################################################
# 8. Function for Spearman correlation
############################################################

calculate_spearman <- function(
    data,
    x,
    y
) {

  tmp <- data[
    !is.na(get(x)) &
      !is.na(get(y))
  ]


  if (nrow(tmp) < 3) {

    return(
      data.table(
        Variable = y,
        N = nrow(tmp),
        Spearman_rho = NA_real_,
        P_value = NA_real_
      )
    )

  }


  test <- cor.test(
    tmp[[x]],
    tmp[[y]],
    method = "spearman",
    exact = FALSE
  )


  data.table(
    Variable = y,
    N = nrow(tmp),
    Spearman_rho = unname(test$estimate),
    P_value = test$p.value
  )

}


############################################################
# 9. Haplotype-level variables
############################################################

hap_variables <- c(
  "INS",
  "DEL",
  "INStelo",
  "DELtelo",
  "INScenSat",
  "DELcenSat"
)


############################################################
# 10. Haplotype-level statistics
############################################################

hap_stats <- rbindlist(
  lapply(
    hap_variables,
    function(v) {

      calculate_spearman(
        hap_df,
        "N50",
        v
      )

    }
  )
)


############################################################
# 11. Add analysis level
############################################################

hap_stats[
  ,
  Level := "Haplotype"
]


############################################################
# 12. Save statistics
############################################################

fwrite(
  hap_stats,
  file.path(
    outdir,
    paste0(
      output_prefix,
      "_haplotype_Spearman.tsv"
    )
  ),
  sep = "\t"
)


############################################################
# 13. Plotting function
############################################################

plot_correlation <- function(
    data,
    x,
    y,
    xlab = "Assembly contig N50",
    ylab = NULL,
    title = NULL
) {

  tmp <- data[
    !is.na(get(x)) &
      !is.na(get(y))
  ]


  test <- cor.test(
    tmp[[x]],
    tmp[[y]],
    method = "spearman",
    exact = FALSE
  )


  rho <- unname(
    test$estimate
  )

  pval <- test$p.value


  if (is.null(ylab)) {
    ylab <- y
  }


  p_label <- ifelse(
    pval < 0.001,
    "P < 0.001",
    paste0(
      "P = ",
      formatC(
        pval,
        format = "e",
        digits = 2
      )
    )
  )


  annotation <- paste0(
    "Spearman \u03c1 = ",
    sprintf("%.3f", rho),
    "\n",
    p_label,
    "\n",
    "n = ",
    nrow(tmp)
  )


  p <- ggplot(
    tmp,
    aes(
      x = .data[[x]],
      y = .data[[y]]
    )
  ) +

    geom_point(
      size = 2.5,
      alpha = 0.7
    ) +

    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.8
    ) +

    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = annotation,
      hjust = 1.05,
      vjust = 1.2,
      size = 4.5
    ) +

    scale_x_continuous(
      labels = label_number(
        scale = 1e-6,
        suffix = " Mb"
      )
    ) +

    labs(
      x = xlab,
      y = ylab,
      title = title
    ) +

    theme_classic(
      base_size = 15
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),

      axis.title = element_text(
        size = 15
      ),

      axis.text = element_text(
        size = 12
      )
    )


  return(p)

}


############################################################
# 14. Generate haplotype-level plots
############################################################

hap_plots <- list()


for (v in hap_variables) {

  hap_plots[[v]] <- plot_correlation(
    hap_df,
    "N50",
    v,
    xlab = "Assembly contig N50",
    ylab = paste0(
      v,
      " variants per haplotype"
    ),
    title = v
  )


  ggsave(
    filename = file.path(
      outdir,
      paste0(
        output_prefix,
        "_N50_vs_",
        v,
        ".pdf"
      )
    ),
    plot = hap_plots[[v]],
    width = 6.5,
    height = 5.5
  )

}


############################################################
# 15. Combine haplotype-level plots
############################################################

hap_combined <- (
  hap_plots$INS +
    hap_plots$DEL
) /

  (
    hap_plots$INStelo +
      hap_plots$DELtelo
  ) /

  (
    hap_plots$INScenSat +
      hap_plots$DELcenSat
  )


ggsave(
  filename = file.path(
    outdir,
    paste0(
      output_prefix,
      "_haplotype_all.pdf"
    )
  ),
  plot = hap_combined,
  width = 13,
  height = 15
)


############################################################
# 16. Sample-level analysis
############################################################
#
# SNP and Indel are duplicated across the two haplotypes.
#
# Therefore they MUST NOT be treated as 260 independent
# observations.
#
# We collapse the assembly N50 information to sample level.
#
# Here we use:
#
# mean N50 of the two haplotypes
#
# as the sample-level N50.
#
############################################################


sample_df <- dat[
  ,
  .(
    N50_mean = mean(
      N50,
      na.rm = TRUE
    ),

    N50_median = median(
      N50,
      na.rm = TRUE
    ),

    N50_min = min(
      N50,
      na.rm = TRUE
    ),

    N50_max = max(
      N50,
      na.rm = TRUE
    ),

    AssemblySize_mean = mean(
      AssemblySize,
      na.rm = TRUE
    ),

    SNP = first(SNP),

    Indel = first(Indel),

    Gender = first(Gender)
  ),
  by = Sample
]


############################################################
# 17. Check duplicated SNP / Indel values
############################################################

dup_check <- dat[
  ,
  .(
    SNP_n = uniqueN(SNP),
    Indel_n = uniqueN(Indel)
  ),
  by = Sample
]


if (any(dup_check$SNP_n > 1)) {

  warning(
    "Some samples have different SNP counts between haplotypes."
  )

}


if (any(dup_check$Indel_n > 1)) {

  warning(
    "Some samples have different Indel counts between haplotypes."
  )

}


############################################################
# 18. Sample-level correlation function
############################################################

sample_variables <- c(
  "SNP",
  "Indel"
)


sample_stats <- rbindlist(
  lapply(
    sample_variables,
    function(v) {

      calculate_spearman(
        sample_df,
        "N50_mean",
        v
      )

    }
  )
)


sample_stats[
  ,
  Level := "Sample"
]


############################################################
# 19. Save sample-level statistics
############################################################

fwrite(
  sample_stats,
  file.path(
    outdir,
    paste0(
      output_prefix,
      "_sample_Spearman.tsv"
    )
  ),
  sep = "\t"
)


############################################################
# 20. Plot SNP vs N50
############################################################

p_snp <- plot_correlation(
  sample_df,
  "N50_mean",
  "SNP",
  xlab = "Mean assembly contig N50 per sample",
  ylab = "SNP discovery count",
  title = "SNP"
)


############################################################
# 21. Plot Indel vs N50
############################################################

p_indel <- plot_correlation(
  sample_df,
  "N50_mean",
  "Indel",
  xlab = "Mean assembly contig N50 per sample",
  ylab = "Indel discovery count",
  title = "Indel"
)


############################################################
# 22. Save sample-level plots
############################################################

ggsave(
  file.path(
    outdir,
    paste0(
      output_prefix,
      "_N50_vs_SNP.pdf"
    )
  ),
  p_snp,
  width = 6.5,
  height = 5.5
)


ggsave(
  file.path(
    outdir,
    paste0(
      output_prefix,
      "_N50_vs_Indel.pdf"
    )
  ),
  p_indel,
  width = 6.5,
  height = 5.5
)


############################################################
# 23. Combined SNP / Indel figure
############################################################

sample_combined <- (
  p_snp +
    p_indel
)


ggsave(
  file.path(
    outdir,
    paste0(
      output_prefix,
      "_sample_SNP_Indel.pdf"
    )
  ),
  sample_combined,
  width = 13,
  height = 5.5
)


############################################################
# 24. Combined statistics
############################################################

all_stats <- rbind(
  hap_stats,
  sample_stats,
  fill = TRUE
)


fwrite(
  all_stats,
  file.path(
    outdir,
    paste0(
      output_prefix,
      "_N50_variant_correlations.tsv"
    )
  ),
  sep = "\t"
)


############################################################
# 25. Print summary
############################################################

cat("\n")
cat("============================================\n")
cat("N50 vs Variant Discovery Analysis\n")
cat("============================================\n\n")


cat("Haplotype-level analysis:\n")
print(hap_stats)


cat("\nSample-level analysis:\n")
print(sample_stats)


cat("\nOutput directory:\n")
cat(normalizePath(outdir), "\n")


############################################################
# END
############################################################

