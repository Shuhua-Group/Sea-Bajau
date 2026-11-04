############################################################
# Assembly quality vs variant discovery
#
# Purpose:
#
# 1. Test whether assembly size is associated with
# downstream variant discovery.
#
# 2. Test whether contig N50 is associated with
# downstream variant discovery.
#
# 3. Evaluate both genome-wide and difficult-to-assemble
# regions (telomeres and centromeres).
#
#
# Analysis levels:
#
# Sample level:
# SNP
# Indel
#
# Haplotype level:
# INS
# DEL
# INStelo
# DELtelo
# INScenSat
# DELcenSat
#
#
# Statistical test:
# Spearman correlation
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

setwd("/path/to/asm_qv")
# Output directory
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
# 5. Check sex composition
############################################################

cat("\n========================================\n")
cat("Sex composition\n")
cat("========================================\n\n")

sex_summary <- df[
  ,
  .(
    N = uniqueN(Sample)
  ),
  by = Gender
]

print(sex_summary)


############################################################
# 6. Create sample-level dataset
############################################################
#
# Important:
#
# SNP and Indel are identical for the two haplotypes
# of the same sample.
#
# Therefore these two variant classes must NOT be
# treated as 260 independent observations.
#
# We retain one row per sample.
#
############################################################

sample_df <- unique(
  df[
    ,
    .(
      Sample,
      Gender,
      N50,
      AssemblySize,
      SNP,
      Indel
    )
  ]
)


############################################################
# 7. Check whether sample-level SNP/Indel values are
# identical between the two haplotypes
############################################################

sample_check <- df[
  ,
  .(
    SNP_unique = uniqueN(SNP),
    Indel_unique = uniqueN(Indel)
  ),
  by = Sample
]


if (
  any(sample_check$SNP_unique > 1) |
  any(sample_check$Indel_unique > 1)
) {

  warning(
    "Some samples have different SNP/Indel values between haplotypes."
  )

}


############################################################
# 8. Haplotype-level dataset
############################################################

hap_df <- df[
  ,
  .(
    Haplotype,
    Sample,
    Gender,
    N50,
    AssemblySize,
    INS,
    DEL,
    INStelo,
    DELtelo,
    INScenSat,
    DELcenSat
  )
]


############################################################
# 9. Variant definitions
############################################################

sample_variants <- c(
  "SNP",
  "Indel"
)


hap_variants <- c(
  "INS",
  "DEL",
  "INStelo",
  "DELtelo",
  "INScenSat",
  "DELcenSat"
)


############################################################
# 10. Function for Spearman correlation
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
# 11. Generate complete statistical table
############################################################

results <- list()

counter <- 1


############################################################
# 11A. AssemblySize
############################################################

for (variant in sample_variants) {

  stat <- calculate_spearman(
    sample_df,
    "AssemblySize",
    variant
  )


  stat[
    ,
    `:=`(
      Metric = "AssemblySize",
      Variant = variant,
      Level = "Sample"
    )
  ]


  results[[counter]] <- stat

  counter <- counter + 1
}


for (variant in hap_variants) {

  stat <- calculate_spearman(
    hap_df,
    "AssemblySize",
    variant
  )


  stat[
    ,
    `:=`(
      Metric = "AssemblySize",
      Variant = variant,
      Level = "Haplotype"
    )
  ]


  results[[counter]] <- stat

  counter <- counter + 1
}


############################################################
# 11B. N50
############################################################

for (variant in sample_variants) {

  stat <- calculate_spearman(
    sample_df,
    "N50",
    variant
  )


  stat[
    ,
    `:=`(
      Metric = "N50",
      Variant = variant,
      Level = "Sample"
    )
  ]


  results[[counter]] <- stat

  counter <- counter + 1
}


for (variant in hap_variants) {

  stat <- calculate_spearman(
    hap_df,
    "N50",
    variant
  )


  stat[
    ,
    `:=`(
      Metric = "N50",
      Variant = variant,
      Level = "Haplotype"
    )
  ]


  results[[counter]] <- stat

  counter <- counter + 1
}


############################################################
# 12. Combine statistical results
############################################################

results_dt <- rbindlist(
  results,
  fill = TRUE
)


setcolorder(
  results_dt,
  c(
    "Metric",
    "Variant",
    "Level",
    "N",
    "rho",
    "P"
  )
)


############################################################
# 13. Multiple-testing correction
############################################################
#
# We have 16 correlation tests:
#
# AssemblySize:
# 8
#
# N50:
# 8
#
# FDR correction is applied across all tests.
#
############################################################

results_dt[
  ,
  FDR := p.adjust(
    P,
    method = "BH"
  )
]


############################################################
# 14. Significance labels
############################################################

results_dt[
  ,
  significance := fifelse(
    P < 0.001,
    "***",
    fifelse(
      P < 0.01,
      "**",
      fifelse(
        P < 0.05,
        "*",
        "ns"
      )
    )
  )
]


############################################################
# 15. Save statistical results
############################################################

fwrite(
  results_dt,
  file.path(
    output_dir,
    "assembly_quality_variant_correlation.tsv"
  ),
  sep = "\t"
)


############################################################
# 16. Print results
############################################################

cat("\n========================================\n")
cat("Correlation results\n")
cat("========================================\n\n")

print(results_dt)


############################################################
# 17. Function for plotting
############################################################

plot_metric <- function(
    data,
    metric,
    variants,
    level
) {


  ########################################################
  # Select correct dataset
  ########################################################

  if (level == "Sample") {

    plot_data <- sample_df[
      ,
      c(
        "Sample",
        metric,
        variants
      ),
      with = FALSE
    ]

  } else {

    plot_data <- hap_df[
      ,
      c(
        "Haplotype",
        "Sample",
        metric,
        variants
      ),
      with = FALSE
    ]

  }


  ########################################################
  # Convert wide → long
  ########################################################

  id_columns <- if (
    level == "Sample"
  ) {

    c(
      "Sample",
      metric
    )

  } else {

    c(
      "Haplotype",
      "Sample",
      metric
    )

  }


  long <- melt(
    plot_data,
    id.vars = id_columns,
    variable.name = "Variant",
    value.name = "VariantCount"
  )


  ########################################################
  # Remove missing values
  ########################################################

  long <- long[
    !is.na(get(metric)) &
      !is.na(VariantCount)
  ]


  ########################################################
  # Variant order
  ########################################################

  long[
    ,
    Variant := factor(
      Variant,
      levels = variants
    )
  ]


  ########################################################
  # Add correlation statistics
  ########################################################

  stat <- results_dt[
    Metric == metric &
      Level == level
  ]


  stat[
    ,
    label := paste0(
      "\u03c1 = ",
      sprintf("%.3f", rho),
      "\nP = ",
      format.pval(
        P,
        digits = 2,
        eps = 1e-3
      ),
      "\nFDR = ",
      format.pval(
        FDR,
        digits = 2,
        eps = 1e-3
      )
    )
  ]


  ########################################################
  # Plot
  ########################################################

  p <- ggplot(
    long,
    aes(
      x = .data[[metric]],
      y = VariantCount
    )
  ) +

    geom_point(
      size = 1.6,
      alpha = 0.55
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
      y = "Variant discovery count"
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


  ########################################################
  # Add correlation labels
  ########################################################

  # Determine annotation positions independently
  # for each facet.

  annotation <- merge(
    long[
      ,
      .(
        xmax = max(
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
    xpos := xmax -
      0.05 *
      (xmax -
         min(
           long[[metric]],
           na.rm = TRUE
         ))
  ]


  annotation[
    ,
    ypos := ymax
  ]


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
# 18. Plot AssemblySize
############################################################

p_assembly_sample <- plot_metric(
  sample_df,
  "AssemblySize",
  sample_variants,
  "Sample"
)


p_assembly_hap <- plot_metric(
  hap_df,
  "AssemblySize",
  hap_variants,
  "Haplotype"
)


############################################################
# 19. Plot N50
############################################################

p_n50_sample <- plot_metric(
  sample_df,
  "N50",
  sample_variants,
  "Sample"
)


p_n50_hap <- plot_metric(
  hap_df,
  "N50",
  hap_variants,
  "Haplotype"
)


############################################################
# 20. Save PDFs
############################################################

ggsave(
  file.path(
    output_dir,
    "AssemblySize_sample_level.pdf"
  ),
  p_assembly_sample,
  width = 12,
  height = 6,
  device = "pdf"
)


ggsave(
  file.path(
    output_dir,
    "AssemblySize_haplotype_level.pdf"
  ),
  p_assembly_hap,
  width = 10,
  height = 12,
  device = "pdf"
)


ggsave(
  file.path(
    output_dir,
    "N50_sample_level.pdf"
  ),
  p_n50_sample,
  width = 12,
  height = 6,
  device = "pdf"
)


ggsave(
  file.path(
    output_dir,
    "N50_haplotype_level.pdf"
  ),
  p_n50_hap,
  width = 10,
  height = 12,
  device = "pdf"
)

############################################################
# 21. End
############################################################

cat("\n========================================\n")
cat("Analysis completed.\n")
cat("Output directory:\n")
cat(output_dir, "\n")
cat("========================================\n")

