############################################################
# Reviewer 2 - Assembly quality and variant discovery
#
# Panel B:
#
# 1. AssemblySize vs variant discovery count
# 2. N50 vs variant discovery count
# 3. AssemblySize vs variant discovery density
# 4. N50 vs variant discovery density
#
# Variant classes:
# SNP
# Indel
# INS
# DEL
#
# IMPORTANT:
# SNP / Indel = sample-level
# INS / DEL = haplotype-level
#
# Correlation:
# Spearman rank correlation
#
# Output:
# - PDF figures
# - correlation statistics
# - analysis table
############################################################


############################################################
# 0. Packages
############################################################

library(data.table)
library(ggplot2)
library(scales)


############################################################
# 1. Parameters
############################################################

# Input TXT file
input_file <- "/path/to/130APC_n50_genomesize.txt"

setwd("/path/to/asm_qv")
# Output directory
output_dir <- "Reviewer2_B_assembly_quality"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 2. Read input
############################################################

dat <- fread(
  input_file,
  sep = "\t",
  header = TRUE
)


############################################################
# 3. Check input columns
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
  "DEL"
)


missing_columns <- setdiff(
  required_columns,
  names(dat)
)


if (length(missing_columns) > 0) {

  stop(
    paste(
      "Missing required columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )

}


############################################################
# 4. Basic data cleaning
############################################################

# Make sure numerical columns are numeric

numeric_columns <- c(
  "N50",
  "AssemblySize",
  "SNP",
  "Indel",
  "INS",
  "DEL"
)


for (x in numeric_columns) {

  dat[
    ,
    (x) := as.numeric(get(x))
  ]

}


# Remove rows with missing assembly metrics

dat <- dat[
  !is.na(N50) &
    !is.na(AssemblySize)
]


############################################################
# 5. Check haplotype structure
############################################################

cat(
  "Number of haplotypes:",
  uniqueN(dat$Haplotype),
  "\n"
)

cat(
  "Number of samples:",
  uniqueN(dat$Sample),
  "\n"
)


############################################################
# 6. Construct sample-level dataset
#
# SNP and Indel are sample-level values.
#
# Therefore, do NOT simply take one haplotype row.
#
# For assembly metrics:
#
# AssemblySize = hap1 + hap2
#
# N50 = mean(hap1, hap2)
#
############################################################

sample_dat <- dat[
  ,
  .(
    N50 = mean(
      N50,
      na.rm = TRUE
    ),

    AssemblySize = sum(
      AssemblySize,
      na.rm = TRUE
    ),

    SNP = first(SNP),

    Indel = first(Indel),

    Gender = first(Gender),

    n_haplotypes = .N
  ),
  by = Sample
]


############################################################
# 7. Check that SNP and Indel are indeed identical
# between the two haplotypes
############################################################

snp_check <- dat[
  ,
  .(
    SNP_n = uniqueN(SNP),
    Indel_n = uniqueN(Indel)
  ),
  by = Sample
]


if (any(snp_check$SNP_n > 1)) {

  warning(
    "Some samples have different SNP values between haplotypes."
  )

}


if (any(snp_check$Indel_n > 1)) {

  warning(
    "Some samples have different Indel values between haplotypes."
  )

}


############################################################
# 8. Calculate sample-level variant density
#
# Density is defined as:
#
# variants / Gb assembled sequence
#
############################################################

sample_dat[
  ,
  SNP_density := SNP /
    (AssemblySize / 1e9)
]


sample_dat[
  ,
  Indel_density := Indel /
    (AssemblySize / 1e9)
]


############################################################
# 9. Construct haplotype-level dataset
#
# INS and DEL are haplotype-specific.
############################################################

hap_dat <- copy(dat)


############################################################
# 10. Calculate haplotype-level SV density
############################################################

hap_dat[
  ,
  INS_density := INS /
    (AssemblySize / 1e9)
]


hap_dat[
  ,
  DEL_density := DEL /
    (AssemblySize / 1e9)
]


############################################################
# 11. Define helper function for Spearman correlation
############################################################

run_spearman <- function(
    data,
    x,
    y,
    metric,
    variant,
    level
) {

  tmp <- data[
    !is.na(get(x)) &
      !is.na(get(y))
  ]


  if (nrow(tmp) < 3) {

    return(
      data.table(
        Level = level,
        Metric = metric,
        Variant = variant,
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
    Level = level,
    Metric = metric,
    Variant = variant,
    N = nrow(tmp),
    Spearman_rho = unname(
      test$estimate
    ),
    P_value = test$p.value
  )

}


############################################################
# 12. Correlation analysis
############################################################

results <- rbindlist(

  list(

    ####################################################
    # Sample level
    #
    # SNP count
    ####################################################

    run_spearman(
      sample_dat,
      "AssemblySize",
      "SNP",
      "AssemblySize",
      "SNP_count",
      "Sample"
    ),

    run_spearman(
      sample_dat,
      "N50",
      "SNP",
      "N50",
      "SNP_count",
      "Sample"
    ),


    ####################################################
    # Sample level
    #
    # Indel count
    ####################################################

    run_spearman(
      sample_dat,
      "AssemblySize",
      "Indel",
      "Indel_count",
      "AssemblySize",
      "Sample"
    ),

    run_spearman(
      sample_dat,
      "N50",
      "Indel",
      "N50",
      "Indel_count",
      "Sample"
    ),


    ####################################################
    # Sample level
    #
    # SNP density
    ####################################################

    run_spearman(
      sample_dat,
      "AssemblySize",
      "SNP_density",
      "AssemblySize",
      "SNP_density",
      "Sample"
    ),

    run_spearman(
      sample_dat,
      "N50",
      "SNP_density",
      "N50",
      "SNP_density",
      "Sample"
    ),


    ####################################################
    # Sample level
    #
    # Indel density
    ####################################################

    run_spearman(
      sample_dat,
      "AssemblySize",
      "Indel_density",
      "AssemblySize",
      "Indel_density",
      "Sample"
    ),

    run_spearman(
      sample_dat,
      "N50",
      "Indel_density",
      "N50",
      "Indel_density",
      "Sample"
    ),


    ####################################################
    # Haplotype level
    #
    # INS count
    ####################################################

    run_spearman(
      hap_dat,
      "AssemblySize",
      "INS",
      "AssemblySize",
      "INS_count",
      "Haplotype"
    ),

    run_spearman(
      hap_dat,
      "N50",
      "INS",
      "N50",
      "INS_count",
      "Haplotype"
    ),


    ####################################################
    # Haplotype level
    #
    # DEL count
    ####################################################

    run_spearman(
      hap_dat,
      "AssemblySize",
      "DEL",
      "AssemblySize",
      "DEL_count",
      "Haplotype"
    ),

    run_spearman(
      hap_dat,
      "N50",
      "DEL",
      "N50",
      "DEL_count",
      "Haplotype"
    ),


    ####################################################
    # Haplotype level
    #
    # INS density
    ####################################################

    run_spearman(
      hap_dat,
      "AssemblySize",
      "INS_density",
      "AssemblySize",
      "INS_density",
      "Haplotype"
    ),

    run_spearman(
      hap_dat,
      "N50",
      "INS_density",
      "N50",
      "INS_density",
      "Haplotype"
    ),


    ####################################################
    # Haplotype level
    #
    # DEL density
    ####################################################

    run_spearman(
      hap_dat,
      "AssemblySize",
      "DEL_density",
      "AssemblySize",
      "DEL_density",
      "Haplotype"
    ),

    run_spearman(
      hap_dat,
      "N50",
      "DEL_density",
      "N50",
      "DEL_density",
      "Haplotype"
    )

  ),
  fill = TRUE
)


############################################################
# 13. Multiple-testing correction
############################################################

results[
  ,
  FDR := p.adjust(
    P_value,
    method = "BH"
  )
]


############################################################
# 14. Save statistics
############################################################

fwrite(
  results,
  file.path(
    output_dir,
    "assembly_quality_variant_discovery_spearman.tsv"
  ),
  sep = "\t"
)


############################################################
# 15. Save processed datasets
############################################################

fwrite(
  sample_dat,
  file.path(
    output_dir,
    "sample_level_analysis.tsv"
  ),
  sep = "\t"
)


fwrite(
  hap_dat,
  file.path(
    output_dir,
    "haplotype_level_analysis.tsv"
  ),
  sep = "\t"
)


############################################################
# 16. Function to generate annotation text
############################################################

get_cor_label <- function(
    data,
    x,
    y
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


  rho <- round(
    unname(test$estimate),
    2
  )


  p <- test$p.value


  if (p < 0.001) {

    p_text <- "P < 0.001"

  } else {

    p_text <- paste0(
      "P = ",
      format.pval(
        p,
        digits = 2,
        eps = 0.001
      )
    )

  }


  paste0(
    "Spearman \u03c1 = ",
    rho,
    "\n",
    p_text
  )

}


############################################################
# 17. General plotting theme
############################################################

theme_review <- theme_classic(
  base_size = 14
) +

  theme(

    axis.title =
      element_text(
        size = 15
      ),

    axis.text =
      element_text(
        size = 12
      ),

    strip.text =
      element_text(
        size = 15,
        face = "bold"
      ),

    plot.title =
      element_text(
        size = 16,
        face = "bold"
      ),

    legend.position =
      "none"

  )


############################################################
# 18. Function for count plots
############################################################

plot_count <- function(
    data,
    x,
    y,
    x_label,
    y_label,
    title
) {


  label <- get_cor_label(
    data,
    x,
    y
  )


  p <- ggplot(
    data,
    aes(
      x = .data[[x]],
      y = .data[[y]]
    )
  ) +

    geom_point(
      size = 2.5,
      alpha = 0.75
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
      label = label,
      hjust = 1.05,
      vjust = 1.2,
      size = 4
    ) +

    labs(
      x = x_label,
      y = y_label,
      title = title
    ) +

    theme_review


  return(p)

}


############################################################
# 19. Function for density plots
############################################################

plot_density <- function(
    data,
    x,
    y,
    x_label,
    y_label,
    title
) {


  label <- get_cor_label(
    data,
    x,
    y
  )


  p <- ggplot(
    data,
    aes(
      x = .data[[x]],
      y = .data[[y]]
    )
  ) +

    geom_point(
      size = 2.5,
      alpha = 0.75
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
      label = label,
      hjust = 1.05,
      vjust = 1.2,
      size = 4
    ) +

    labs(
      x = x_label,
      y = y_label,
      title = title
    ) +

    theme_review


  return(p)

}


############################################################
# 20. AssemblySize vs SNP count
############################################################

p_assembly_snp_count <- plot_count(

  sample_dat,

  "AssemblySize",
  "SNP",

  "Assembly size (bp)",
  "SNP count",

  "Assembly size vs. SNP discovery"
)


ggsave(
  file.path(
    output_dir,
    "B_AssemblySize_vs_SNP_count.pdf"
  ),
  p_assembly_snp_count,
  width = 6,
  height = 5
)


############################################################
# 21. N50 vs SNP count
############################################################

p_n50_snp_count <- plot_count(

  sample_dat,

  "N50",
  "SNP",

  "Contig N50 (bp)",
  "SNP count",

  "Contig N50 vs. SNP discovery"
)


ggsave(
  file.path(
    output_dir,
    "B_N50_vs_SNP_count.pdf"
  ),
  p_n50_snp_count,
  width = 6,
  height = 5
)


############################################################
# 22. AssemblySize vs Indel count
############################################################

p_assembly_indel_count <- plot_count(

  sample_dat,

  "AssemblySize",
  "Indel",

  "Assembly size (bp)",
  "Indel count",

  "Assembly size vs. indel discovery"
)


ggsave(
  file.path(
    output_dir,
    "B_AssemblySize_vs_Indel_count.pdf"
  ),
  p_assembly_indel_count,
  width = 6,
  height = 5
)


############################################################
# 23. N50 vs Indel count
############################################################

p_n50_indel_count <- plot_count(

  sample_dat,

  "N50",
  "Indel",

  "Contig N50 (bp)",
  "Indel count",

  "Contig N50 vs. indel discovery"
)


ggsave(
  file.path(
    output_dir,
    "B_N50_vs_Indel_count.pdf"
  ),
  p_n50_indel_count,
  width = 6,
  height = 5
)


############################################################
# 24. AssemblySize vs SNP density
############################################################

p_assembly_snp_density <- plot_density(

  sample_dat,

  "AssemblySize",
  "SNP_density",

  "Assembly size (bp)",
  "SNP density (variants / Gb)",

  "Assembly size vs. SNP discovery density"
)


ggsave(
  file.path(
    output_dir,
    "B_AssemblySize_vs_SNP_density.pdf"
  ),
  p_assembly_snp_density,
  width = 6,
  height = 5
)


############################################################
# 25. N50 vs SNP density
############################################################

p_n50_snp_density <- plot_density(

  sample_dat,

  "N50",
  "SNP_density",

  "Contig N50 (bp)",
  "SNP density (variants / Gb)",

  "Contig N50 vs. SNP discovery density"
)


ggsave(
  file.path(
    output_dir,
    "B_N50_vs_SNP_density.pdf"
  ),
  p_n50_snp_density,
  width = 6,
  height = 5
)


############################################################
# 26. AssemblySize vs Indel density
############################################################

p_assembly_indel_density <- plot_density(

  sample_dat,

  "AssemblySize",
  "Indel_density",

  "Assembly size (bp)",
  "Indel density (variants / Gb)",

  "Assembly size vs. indel discovery density"
)


ggsave(
  file.path(
    output_dir,
    "B_AssemblySize_vs_Indel_density.pdf"
  ),
  p_assembly_indel_density,
  width = 6,
  height = 5
)


############################################################
# 27. N50 vs Indel density
############################################################

p_n50_indel_density <- plot_density(

  sample_dat,

  "N50",
  "Indel_density",

  "Contig N50 (bp)",
  "Indel density (variants / Gb)",

  "Contig N50 vs. indel discovery density"
)


ggsave(
  file.path(
    output_dir,
    "B_N50_vs_Indel_density.pdf"
  ),
  p_n50_indel_density,
  width = 6,
  height = 5
)


############################################################
# 28. AssemblySize vs INS count
############################################################

p_assembly_ins_count <- plot_count(

  hap_dat,

  "AssemblySize",
  "INS",

  "Assembly size (bp)",
  "Insertion count",

  "Assembly size vs. insertion discovery"
)


ggsave(
  file.path(
    output_dir,
    "B_AssemblySize_vs_INS_count.pdf"
  ),
  p_assembly_ins_count,
  width = 6,
  height = 5
)


############################################################
# 29. N50 vs INS count
############################################################

p_n50_ins_count <- plot_count(

  hap_dat,

  "N50",
  "INS",

  "Contig N50 (bp)",
  "Insertion count",

  "Contig N50 vs. insertion discovery"
)


ggsave(
  file.path(
    output_dir,
    "B_N50_vs_INS_count.pdf"
  ),
  p_n50_ins_count,
  width = 6,
  height = 5
)


############################################################
# 30. AssemblySize vs DEL count
############################################################

p_assembly_del_count <- plot_count(

  hap_dat,

  "AssemblySize",
  "DEL",

  "Assembly size (bp)",
  "Deletion count",

  "Assembly size vs. deletion discovery"
)


ggsave(
  file.path(
    output_dir,
    "B_AssemblySize_vs_DEL_count.pdf"
  ),
  p_assembly_del_count,
  width = 6,
  height = 5
)


############################################################
# 31. N50 vs DEL count
############################################################

p_n50_del_count <- plot_count(

  hap_dat,

  "N50",
  "DEL",

  "Contig N50 (bp)",
  "Deletion count",

  "Contig N50 vs. deletion discovery"
)


ggsave(
  file.path(
    output_dir,
    "B_N50_vs_DEL_count.pdf"
  ),
  p_n50_del_count,
  width = 6,
  height = 5
)


############################################################
# 32. AssemblySize vs INS density
############################################################

p_assembly_ins_density <- plot_density(

  hap_dat,

  "AssemblySize",
  "INS_density",

  "Assembly size (bp)",
  "Insertion density (variants / Gb)",

  "Assembly size vs. insertion discovery density"
)


ggsave(
  file.path(
    output_dir,
    "B_AssemblySize_vs_INS_density.pdf"
  ),
  p_assembly_ins_density,
  width = 6,
  height = 5
)


############################################################
# 33. N50 vs INS density
############################################################

p_n50_ins_density <- plot_density(

  hap_dat,

  "N50",
  "INS_density",

  "Contig N50 (bp)",
  "Insertion density (variants / Gb)",

  "Contig N50 vs. insertion discovery density"
)


ggsave(
  file.path(
    output_dir,
    "B_N50_vs_INS_density.pdf"
  ),
  p_n50_ins_density,
  width = 6,
  height = 5
)


############################################################
# 34. AssemblySize vs DEL density
############################################################

p_assembly_del_density <- plot_density(

  hap_dat,

  "AssemblySize",
  "DEL_density",

  "Assembly size (bp)",
  "Deletion density (variants / Gb)",

  "Assembly size vs. deletion discovery density"
)


ggsave(
  file.path(
    output_dir,
    "B_AssemblySize_vs_DEL_density.pdf"
  ),
  p_assembly_del_density,
  width = 6,
  height = 5
)


############################################################
# 35. N50 vs DEL density
############################################################

p_n50_del_density <- plot_density(

  hap_dat,

  "N50",
  "DEL_density",

  "Contig N50 (bp)",
  "Deletion density (variants / Gb)",

  "Contig N50 vs. deletion discovery density"
)


ggsave(
  file.path(
    output_dir,
    "B_N50_vs_DEL_density.pdf"
  ),
  p_n50_del_density,
  width = 6,
  height = 5
)


############################################################
# 36. Finish
############################################################

cat("\n")
cat("============================================\n")
cat("Analysis completed.\n")
cat("============================================\n")
cat(
  "Results:",
  output_dir,
  "\n"
)
cat(
  "Correlation table:",
  file.path(
    output_dir,
    "assembly_quality_variant_discovery_spearman.tsv"
  ),
  "\n"
)
cat("============================================\n")

