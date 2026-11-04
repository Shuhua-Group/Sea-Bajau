############################################################
# Assembly quality vs variant discovery
#
# Updated analysis design
#
# Statistical unit:
# Individual / sample level
#
# Number of individuals:
# 130
#
# Each individual has two haplotypes.
#
# Sample-level variables:
# SNP
# Indel
#
# Haplotype-derived variables aggregated to individual level:
# INS
# DEL
# INStelo
# DELtelo
# INScenSat
# DELcenSat
#
# Assembly quality:
# AssemblySize = sum of two haplotypes
# N50 = mean of two haplotype N50 values
#
# Covariate:
# Gender
#
# Primary statistical tests:
# 1. Spearman correlation
# 2. Partial Spearman correlation controlling for Gender
#
# Sensitivity analysis:
# Male only
# Female only
#
############################################################


############################################################
# 0. Load packages
############################################################

library(data.table)
library(ggplot2)
library(scales)


############################################################
# 1. Parameters
############################################################

input_file <- "/path/to/130APC_n50_genomesize.txt"

setwd(
  "/path/to/asm_qv"
)

output_dir <- "assembly_quality_analysis"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


############################################################
# 2. Read data
############################################################

df <- fread(
  input_file,
  header = TRUE,
  sep = "\t"
)


############################################################
# 3. Basic data checking
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
  colnames(df)
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
# 4. Data type conversion
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

for (x in numeric_columns) {

  df[
    ,
    (x) := as.numeric(get(x))
  ]

}


############################################################
# 5. Basic sample / haplotype checks
############################################################

cat("\n========================================\n")
cat("Basic data information\n")
cat("========================================\n\n")

cat(
  "Total rows:",
  nrow(df),
  "\n"
)

cat(
  "Unique samples:",
  uniqueN(df$Sample),
  "\n"
)

cat(
  "Unique haplotypes:",
  uniqueN(df$Haplotype),
  "\n"
)


############################################################
# 6. Check number of haplotypes per sample
############################################################

hap_count_check <- df[
  ,
  .(
    N_haplotype = uniqueN(Haplotype)
  ),
  by = Sample
]

cat("\n========================================\n")
cat("Haplotype count per sample\n")
cat("========================================\n\n")

print(
  hap_count_check[
    ,
    .(
      N_samples = .N
    ),
    by = N_haplotype
  ]
)

if (any(hap_count_check$N_haplotype != 2)) {

  warning(
    "Some samples do not have exactly two haplotypes."
  )

}


############################################################
# 7. Check Gender consistency within sample
############################################################

gender_check <- df[
  ,
  .(
    Gender_unique = uniqueN(Gender)
  ),
  by = Sample
]

if (any(gender_check$Gender_unique > 1)) {

  stop(
    "Some samples have inconsistent Gender annotations."
  )

}


############################################################
# 8. Check sex composition
############################################################

cat("\n========================================\n")
cat("Sex composition\n")
cat("========================================\n\n")

sex_summary <- unique(
  df[
    ,
    .(
      Sample,
      Gender
    )
  ]
)[
  ,
  .(
    N = .N
  ),
  by = Gender
]

print(sex_summary)


############################################################
# 9. Check whether SNP / Indel are identical between
# the two haplotypes
############################################################

sample_variant_check <- df[
  ,
  .(
    N50_unique = uniqueN(N50),
    AssemblySize_unique = uniqueN(AssemblySize),
    SNP_unique = uniqueN(SNP),
    Indel_unique = uniqueN(Indel)
  ),
  by = Sample
]

cat("\n========================================\n")
cat("Within-sample consistency check\n")
cat("========================================\n\n")

cat(
  "Samples with different SNP values between haplotypes:",
  sum(sample_variant_check$SNP_unique > 1),
  "\n"
)

cat(
  "Samples with different Indel values between haplotypes:",
  sum(sample_variant_check$Indel_unique > 1),
  "\n"
)

cat(
  "Samples with different N50 values between haplotypes:",
  sum(sample_variant_check$N50_unique > 1),
  "\n"
)

cat(
  "Samples with different AssemblySize values between haplotypes:",
  sum(sample_variant_check$AssemblySize_unique > 1),
  "\n"
)


if (any(sample_variant_check$SNP_unique > 1)) {

  warning(
    "Some samples have different SNP values between haplotypes."
  )

}

if (any(sample_variant_check$Indel_unique > 1)) {

  warning(
    "Some samples have different Indel values between haplotypes."
  )

}


############################################################
# 10. Create individual-level dataset
############################################################
#
# IMPORTANT
#
# Each individual contributes ONE observation.
#
# SNP / Indel:
# Sample-level values.
#
# AssemblySize:
# Sum of the two haplotype assembly sizes.
#
# N50:
# Mean of the two haplotype N50 values.
#
# SV counts:
# Sum across the two haplotypes.
#
# Gender:
# One value per individual.
#
############################################################

sample_df <- df[
  ,
  .(

    Gender = first(Gender),

    ########################################################
    # Assembly quality
    ########################################################

    AssemblySize = sum(
      AssemblySize,
      na.rm = TRUE
    ),

    N50 = mean(
      N50,
      na.rm = TRUE
    ),

    ########################################################
    # Sample-level small variants
    ########################################################

    SNP = first(SNP),

    Indel = first(Indel),

    ########################################################
    # Diploid / individual-level SV discovery
    ########################################################

    INS = sum(
      INS,
      na.rm = TRUE
    ),

    DEL = sum(
      DEL,
      na.rm = TRUE
    ),

    INStelo = sum(
      INStelo,
      na.rm = TRUE
    ),

    DELtelo = sum(
      DELtelo,
      na.rm = TRUE
    ),

    INScenSat = sum(
      INScenSat,
      na.rm = TRUE
    ),

    DELcenSat = sum(
      DELcenSat,
      na.rm = TRUE
    )

  ),
  by = Sample
]


############################################################
# 11. Verify final sample-level dataset
############################################################

cat("\n========================================\n")
cat("Final sample-level dataset\n")
cat("========================================\n\n")

cat(
  "Number of samples:",
  nrow(sample_df),
  "\n"
)

cat(
  "Number of unique samples:",
  uniqueN(sample_df$Sample),
  "\n"
)

if (nrow(sample_df) != 130) {

  warning(
    paste(
      "Expected 130 samples, but found",
      nrow(sample_df)
    )
  )

}


############################################################
# 12. Verify Gender after aggregation
############################################################

cat("\n========================================\n")
cat("Gender composition after aggregation\n")
cat("========================================\n\n")

print(
  sample_df[
    ,
    .(
      N = .N
    ),
    by = Gender
  ]
)


############################################################
# 13. Define variant groups
############################################################

sample_variants <- c(
  "SNP",
  "Indel",
  "INS",
  "DEL",
  "INStelo",
  "DELtelo",
  "INScenSat",
  "DELcenSat"
)


############################################################
# 14. Function for ordinary Spearman correlation
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
        N = nrow(tmp),
        rho = NA_real_,
        P = NA_real_
      )
    )

  }

  test <- suppressWarnings(
    cor.test(
      tmp[[x]],
      tmp[[y]],
      method = "spearman",
      exact = FALSE
    )
  )

  data.table(
    N = nrow(tmp),
    rho = unname(test$estimate),
    P = test$p.value
  )

}


############################################################
# 15. Function for partial Spearman correlation
#
# Controls for Gender.
#
# Procedure:
#
# 1. Rank-transform X
# 2. Rank-transform Y
# 3. Regress ranked X against Gender
# 4. Regress ranked Y against Gender
# 5. Correlate the two residuals
#
############################################################

calculate_partial_spearman <- function(
    data,
    x,
    y,
    covariate = "Gender"
) {

  tmp <- data[
    !is.na(get(x)) &
      !is.na(get(y)) &
      !is.na(get(covariate))
  ]

  if (nrow(tmp) < 5) {

    return(
      data.table(
        N = nrow(tmp),
        rho_partial = NA_real_,
        P_partial = NA_real_
      )
    )

  }

  ##########################################################
  # Convert Gender to factor
  ##########################################################

  tmp[
    ,
    Gender_factor := factor(
      get(covariate)
    )
  ]

  if (length(unique(tmp$Gender_factor)) < 2) {

    return(
      data.table(
        N = nrow(tmp),
        rho_partial = NA_real_,
        P_partial = NA_real_
      )
    )

  }

  ##########################################################
  # Rank transform
  ##########################################################

  x_rank <- rank(
    tmp[[x]],
    ties.method = "average"
  )

  y_rank <- rank(
    tmp[[y]],
    ties.method = "average"
  )

  ##########################################################
  # Remove Gender effect
  ##########################################################

  x_resid <- residuals(
    lm(
      x_rank ~ Gender_factor,
      data = tmp
    )
  )

  y_resid <- residuals(
    lm(
      y_rank ~ Gender_factor,
      data = tmp
    )
  )

  ##########################################################
  # Correlation of residuals
  ##########################################################

  test <- suppressWarnings(
    cor.test(
      x_resid,
      y_resid,
      method = "pearson"
    )
  )

  data.table(
    N = nrow(tmp),
    rho_partial = unname(test$estimate),
    P_partial = test$p.value
  )

}


############################################################
# 16. Primary analysis
############################################################
#
# 8 variants × 2 assembly metrics = 16 tests
#
# FDR correction is performed across these 16 tests
# separately for:
#
# A. Ordinary Spearman
# B. Gender-adjusted partial Spearman
#
############################################################

primary_results <- list()

counter <- 1


############################################################
# 16A. AssemblySize
############################################################

for (variant in sample_variants) {

  stat <- calculate_spearman(
    sample_df,
    "AssemblySize",
    variant
  )

  partial <- calculate_partial_spearman(
    sample_df,
    "AssemblySize",
    variant,
    "Gender"
  )

  stat[
    ,
    `:=`(
      Metric = "AssemblySize",
      Variant = variant,
      Level = "Sample",
      rho_partial = partial$rho_partial,
      P_partial = partial$P_partial
    )
  ]

  primary_results[[counter]] <- stat

  counter <- counter + 1

}


############################################################
# 16B. N50
############################################################

for (variant in sample_variants) {

  stat <- calculate_spearman(
    sample_df,
    "N50",
    variant
  )

  partial <- calculate_partial_spearman(
    sample_df,
    "N50",
    variant,
    "Gender"
  )

  stat[
    ,
    `:=`(
      Metric = "N50",
      Variant = variant,
      Level = "Sample",
      rho_partial = partial$rho_partial,
      P_partial = partial$P_partial
    )
  ]

  primary_results[[counter]] <- stat

  counter <- counter + 1

}


############################################################
# 17. Combine primary results
############################################################

primary_results_dt <- rbindlist(
  primary_results,
  fill = TRUE
)

setcolorder(
  primary_results_dt,
  c(
    "Metric",
    "Variant",
    "Level",
    "N",
    "rho",
    "P",
    "rho_partial",
    "P_partial"
  )
)


############################################################
# 18. FDR correction
############################################################

primary_results_dt[
  ,
  FDR := p.adjust(
    P,
    method = "BH"
  )
]

primary_results_dt[
  ,
  FDR_partial := p.adjust(
    P_partial,
    method = "BH"
  )
]


############################################################
# 19. Significance labels
############################################################

primary_results_dt[
  ,
  significance := fifelse(
    FDR < 0.001,
    "***",
    fifelse(
      FDR < 0.01,
      "**",
      fifelse(
        FDR < 0.05,
        "*",
        "ns"
      )
    )
  )
]

primary_results_dt[
  ,
  significance_partial := fifelse(
    FDR_partial < 0.001,
    "***",
    fifelse(
      FDR_partial < 0.01,
      "**",
      fifelse(
        FDR_partial < 0.05,
        "*",
        "ns"
      )
    )
  )
]


############################################################
# 20. Save primary results
############################################################

fwrite(
  primary_results_dt,
  file.path(
    output_dir,
    "assembly_quality_variant_correlation_gender_adjusted.tsv"
  ),
  sep = "\t"
)


############################################################
# 21. Print primary results
############################################################

cat("\n========================================\n")
cat("Primary correlation results\n")
cat("========================================\n\n")

print(
  primary_results_dt
)


############################################################
# 22. Male-only and Female-only sensitivity analysis
############################################################

calculate_sex_specific <- function(
    data,
    sex
) {

  tmp <- data[
    Gender == sex
  ]

  result_list <- list()

  counter <- 1

  for (metric in c(
    "AssemblySize",
    "N50"
  )) {

    for (variant in sample_variants) {

      stat <- calculate_spearman(
        tmp,
        metric,
        variant
      )

      stat[
        ,
        `:=`(
          Gender = sex,
          Metric = metric,
          Variant = variant,
          Level = "Sample"
        )
      ]

      result_list[[counter]] <- stat

      counter <- counter + 1

    }

  }

  rbindlist(
    result_list,
    fill = TRUE
  )

}


############################################################
# 23. Run sex-specific analyses
############################################################

male_results <- calculate_sex_specific(
  sample_df,
  "Male"
)

female_results <- calculate_sex_specific(
  sample_df,
  "Female"
)

sex_specific_results <- rbindlist(
  list(
    male_results,
    female_results
  ),
  fill = TRUE
)


############################################################
# 24. Save sex-specific results
############################################################

fwrite(
  sex_specific_results,
  file.path(
    output_dir,
    "assembly_quality_variant_correlation_sex_specific.tsv"
  ),
  sep = "\t"
)


############################################################
# 25. Print sex-specific results
############################################################

cat("\n========================================\n")
cat("Sex-specific correlation results\n")
cat("========================================\n\n")

print(
  sex_specific_results
)


############################################################
# 26. Plotting function
############################################################

plot_metric <- function(
    data,
    metric,
    variants,
    adjusted = TRUE
) {

  ##########################################################
  # Prepare data
  ##########################################################

  plot_data <- data[
    ,
    c(
      "Sample",
      "Gender",
      metric,
      variants
    ),
    with = FALSE
  ]


  ##########################################################
  # Wide → long
  ##########################################################

  long <- melt(
    plot_data,
    id.vars = c(
      "Sample",
      "Gender",
      metric
    ),
    variable.name = "Variant",
    value.name = "VariantCount"
  )


  ##########################################################
  # Remove missing values
  ##########################################################

  long <- long[
    !is.na(get(metric)) &
      !is.na(VariantCount)
  ]


  ##########################################################
  # Variant order
  ##########################################################

  long[
    ,
    Variant := factor(
      Variant,
      levels = variants
    )
  ]


  ##########################################################
  # Correlation statistics
  ##########################################################

  stat <- primary_results_dt[
    Metric == metric
  ]


  ##########################################################
  # Labels
  ##########################################################

  if (adjusted) {

    stat[
      ,
      label := paste0(
        "\u03c1 = ",
        sprintf(
          "%.3f",
          rho_partial
        ),
        "\nFDR = ",
        format.pval(
          FDR_partial,
          digits = 2,
          eps = 1e-3
        )
      )
    ]

  } else {

    stat[
      ,
      label := paste0(
        "\u03c1 = ",
        sprintf(
          "%.3f",
          rho
        ),
        "\nFDR = ",
        format.pval(
          FDR,
          digits = 2,
          eps = 1e-3
        )
      )
    ]

  }


  ##########################################################
  # Main plot
  ##########################################################

  p <- ggplot(
    long,
    aes(
      x = .data[[metric]],
      y = VariantCount
    )
  ) +

    geom_point(
      aes(
        shape = Gender
      ),
      size = 1.8,
      alpha = 0.60
    ) +

    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.7
    ) +

    facet_wrap(
      ~Variant,
      scales = "free_y",
      ncol = 2
    ) +

    scale_x_continuous(
      labels = label_number(
        scale = 1e-9,
        suffix = " Gb"
      )
    ) +

    labs(
      x = metric,
      y = "Variant discovery count",
      shape = "Gender"
    ) +

    theme_classic(
      base_size = 15
    ) +

    theme(
      strip.text = element_text(
        face = "bold",
        size = 14
      ),

      axis.title = element_text(
        face = "bold"
      ),

      panel.spacing = unit(
        1.2,
        "lines"
      )
    )


  ##########################################################
  # Annotation positions
  ##########################################################

  annotation <- merge(
    long[
      ,
      .(
        xmax = max(
          get(metric),
          na.rm = TRUE
        ),
        xmin = min(
          get(metric),
          na.rm = TRUE
        ),
        ymax = max(
          VariantCount,
          na.rm = TRUE
        )
      ),
      by = Variant
    ],
    stat[
      ,
      .(
        Variant,
        label
      )
    ],
    by = "Variant",
    all.x = TRUE
  )


  annotation[
    ,
    xpos :=
      xmax -
      0.05 *
      (xmax - xmin)
  ]


  annotation[
    ,
    ypos := ymax
  ]


  ##########################################################
  # Add labels
  ##########################################################

  p <- p +

    geom_text(
      data = annotation,
      aes(
        x = xpos,
        y = ypos,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.5
    )


  return(p)

}


############################################################
# 27. Generate plots
############################################################

p_assembly <- plot_metric(
  sample_df,
  "AssemblySize",
  sample_variants,
  adjusted = TRUE
)

p_n50 <- plot_metric(
  sample_df,
  "N50",
  sample_variants,
  adjusted = TRUE
)


############################################################
# 28. Save plots
############################################################

ggsave(
  file.path(
    output_dir,
    "AssemblySize_sample_level_gender_adjusted.pdf"
  ),
  p_assembly,
  width = 12,
  height = 12,
  device = "pdf"
)

ggsave(
  file.path(
    output_dir,
    "N50_sample_level_gender_adjusted.pdf"
  ),
  p_n50,
  width = 12,
  height = 12,
  device = "pdf"
)


############################################################
# 29. Save final sample-level dataset
############################################################

fwrite(
  sample_df,
  file.path(
    output_dir,
    "sample_level_assembly_variant_data.tsv"
  ),
  sep = "\t"
)


############################################################
# 30. Final summary
############################################################

cat("\n========================================\n")
cat("Analysis completed.\n")
cat("========================================\n\n")

cat(
  "Number of individuals:",
  nrow(sample_df),
  "\n"
)

cat(
  "Number of haplotypes in input:",
  uniqueN(df$Haplotype),
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)

cat("\nFiles generated:\n")

cat(
  "1. sample_level_assembly_variant_data.tsv\n"
)

cat(
  "2. assembly_quality_variant_correlation_gender_adjusted.tsv\n"
)

cat(
  "3. assembly_quality_variant_correlation_sex_specific.tsv\n"
)

cat(
  "4. AssemblySize_sample_level_gender_adjusted.pdf\n"
)

cat(
  "5. N50_sample_level_gender_adjusted.pdf\n"
)

cat("\n========================================\n")

