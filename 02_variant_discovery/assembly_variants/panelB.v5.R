############################################################
# Assembly quality vs variant discovery
#
# Revised analysis
#
# Main question:
#
# Does assembly quality affect SV discovery after accounting
# for population differences?
#
# Analysis strategy:
#
# 1. Collapse two haplotypes into one sample
# 2. Merge haplotype-level SV counts
# 3. Keep Male and Female separate
# 4. Descriptive Spearman correlation
# 5. Population-adjusted negative binomial regression
#
# Models:
#
# SV_count ~ AssemblySize + Population
# SV_count ~ N50 + Population
#
# No variant density is used.
#
############################################################


############################################################
# 0. Load packages
############################################################

library(data.table)
library(ggplot2)
library(scales)
library(MASS)


############################################################
# 1. Parameters
############################################################

input_file <- "/path/to/130APC_n50_genomesize.txt"

setwd(
  "/path/to/asm_qv"
)


output_dir <- "assembly_quality_population_adjusted"

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
# 3. Required columns
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
  "DELcenSat",
  "Population"
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
# 4. Convert numeric columns
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
# 5. Basic checks
############################################################

cat("\n========================================\n")
cat("Basic information\n")
cat("========================================\n\n")

cat("Total rows:", nrow(df), "\n")
cat("Unique samples:", uniqueN(df$Sample), "\n")
cat("Unique haplotypes:", uniqueN(df$Haplotype), "\n\n")


cat("Gender composition:\n")
print(
  df[
    ,
    .(
      N_samples = uniqueN(Sample)
    ),
    by = Gender
  ]
)


cat("\nPopulation composition:\n")
print(
  df[
    ,
    .(
      N_samples = uniqueN(Sample)
    ),
    by = Population
  ][order(-N_samples)]
)


############################################################
# 6. Check that each sample has two haplotypes
############################################################

hap_check <- df[
  ,
  .(
    N_haplotypes = uniqueN(Haplotype),
    N_gender = uniqueN(Gender),
    N_population = uniqueN(Population)
  ),
  by = Sample
]


if (any(hap_check$N_haplotypes != 2)) {

  warning(
    "Some samples do not have exactly two haplotypes."
  )

}


if (any(hap_check$N_gender != 1)) {

  stop(
    "Some samples have inconsistent Gender assignments."
  )
}


if (any(hap_check$N_population != 1)) {

  stop(
    "Some samples have inconsistent Population assignments."
  )
}


############################################################
# 7. Check SNP / Indel consistency
############################################################
#
# SNP and Indel are sample-level quantities.
#
# They should be identical for the two haplotypes.
#
############################################################

sample_variant_check <- df[
  ,
  .(
    SNP_unique = uniqueN(SNP),
    Indel_unique = uniqueN(Indel)
  ),
  by = Sample
]


if (
  any(sample_variant_check$SNP_unique > 1) |
  any(sample_variant_check$Indel_unique > 1)
) {

  warning(
    "Some samples have different SNP/Indel values between haplotypes."
  )
}


############################################################
# 8. Collapse haplotypes into samples
############################################################
#
# IMPORTANT
#
# All SV counts are summed across the two haplotypes.
#
# Assembly quality:
#
# AssemblySize = hap1 + hap2
#
# N50 = mean(hap1, hap2)
#
# SNP / Indel:
#
# retain one copy because these are sample-level values.
#
############################################################

sample_df <- df[
  ,
  .(

    Gender = first(Gender),

    Population = first(Population),

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
    # Sample-level variants
    ########################################################

    SNP = first(SNP),

    Indel = first(Indel),

    ########################################################
    # Haplotype-level SVs → sample-level SV burden
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
# 9. Convert categorical variables to factors
############################################################

sample_df[
  ,
  Gender := factor(Gender)
]


sample_df[
  ,
  Population := factor(Population)
]


############################################################
# 10. Check final sample numbers
############################################################

cat("\n========================================\n")
cat("Collapsed sample-level dataset\n")
cat("========================================\n\n")

cat(
  "Total samples:",
  nrow(sample_df),
  "\n"
)


cat("\nSamples by gender:\n")

print(
  sample_df[
    ,
    .(
      N = .N
    ),
    by = Gender
  ]
)


cat("\nSamples by gender and population:\n")

print(
  sample_df[
    ,
    .(
      N = .N
    ),
    by = .(Gender, Population)
  ][order(Gender, -N)]
)


############################################################
# 11. Save collapsed dataset
############################################################

fwrite(
  sample_df,
  file.path(
    output_dir,
    "sample_level_collapsed_dataset.tsv"
  ),
  sep = "\t"
)


############################################################
# 12. Variant definitions
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


sv_variants <- c(
  "INS",
  "DEL",
  "INStelo",
  "DELtelo",
  "INScenSat",
  "DELcenSat"
)


############################################################
# 13. Function:
# Spearman correlation
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
# 14. Descriptive Spearman analysis
############################################################
#
# Male and Female are analyzed separately.
#
############################################################

spearman_results <- list()

counter <- 1


for (gender_i in levels(sample_df$Gender)) {

  gender_df <- sample_df[
    Gender == gender_i
  ]


  ##########################################################
  # AssemblySize
  ##########################################################

  for (variant in sample_variants) {

    stat <- calculate_spearman(
      gender_df,
      "AssemblySize",
      variant
    )


    stat[
      ,
      `:=`(
        Gender = gender_i,
        Metric = "AssemblySize",
        Variant = variant,
        Level = "Sample"
      )
    ]


    spearman_results[[counter]] <- stat

    counter <- counter + 1
  }


  ##########################################################
  # N50
  ##########################################################

  for (variant in sample_variants) {

    stat <- calculate_spearman(
      gender_df,
      "N50",
      variant
    )


    stat[
      ,
      `:=`(
        Gender = gender_i,
        Metric = "N50",
        Variant = variant,
        Level = "Sample"
      )
    ]


    spearman_results[[counter]] <- stat

    counter <- counter + 1
  }
}


spearman_dt <- rbindlist(
  spearman_results,
  fill = TRUE
)


setcolorder(
  spearman_dt,
  c(
    "N",
    "rho",
    "P",
    "Gender",
    "Metric",
    "Variant",
    "Level"
  )
)


############################################################
# 15. FDR correction for descriptive analysis
############################################################

spearman_dt[
  ,
  FDR := p.adjust(
    P,
    method = "BH"
  ),
  by = Gender
]


spearman_dt[
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


############################################################
# 16. Save Spearman results
############################################################

fwrite(
  spearman_dt,
  file.path(
    output_dir,
    "01_descriptive_spearman.tsv"
  ),
  sep = "\t"
)


############################################################
# 17. Print Spearman results
############################################################

cat("\n========================================\n")
cat("Descriptive Spearman correlations\n")
cat("========================================\n\n")

print(
  spearman_dt
)


############################################################
# 18. Population-adjusted negative binomial model
############################################################
#
# Main analysis
#
# For each Gender:
#
# Variant ~ AssemblySize + Population
#
# Variant ~ N50 + Population
#
# The coefficient of AssemblySize / N50 tells us whether
# assembly quality explains additional variation in
# variant counts after accounting for Population.
#
############################################################

run_nb_model <- function(
    data,
    variant,
    metric
) {

  tmp <- data[
    !is.na(get(variant)) &
      !is.na(get(metric)) &
      !is.na(Population)
  ]


  ########################################################
  # Check number of populations
  ########################################################

  n_pop <- uniqueN(
    tmp$Population
  )


  if (nrow(tmp) < 10) {

    return(
      data.table(
        N = nrow(tmp),
        N_population = n_pop,
        Estimate = NA_real_,
        SE = NA_real_,
        Z = NA_real_,
        P = NA_real_,
        IRR = NA_real_,
        Model = "NegativeBinomial",
        Status = "Too few samples"
      )
    )

  }


  ########################################################
  # Check whether population has sufficient replication
  ########################################################

  pop_table <- table(
    tmp$Population
  )


  if (any(pop_table < 2)) {

    warning(
      paste(
        "Some populations have only one sample for",
        variant,
        metric
      )
    )

  }


  ########################################################
  # Construct formula
  ########################################################

  formula_text <- paste(
    variant,
    "~",
    metric,
    "+ Population"
  )


  formula <- as.formula(
    formula_text
  )


  ########################################################
  # Negative binomial model
  ########################################################

  fit <- tryCatch(

    glm.nb(
      formula,
      data = tmp
    ),

    error = function(e) {
      NULL
    }

  )


  ########################################################
  # If model fails
  ########################################################

  if (is.null(fit)) {

    return(
      data.table(
        N = nrow(tmp),
        N_population = n_pop,
        Estimate = NA_real_,
        SE = NA_real_,
        Z = NA_real_,
        P = NA_real_,
        IRR = NA_real_,
        Model = "NegativeBinomial",
        Status = "Model failed"
      )
    )

  }


  ########################################################
  # Extract coefficient
  ########################################################

  coef_table <- summary(
    fit
  )$coefficients


  if (!(metric %in% rownames(coef_table))) {

    return(
      data.table(
        N = nrow(tmp),
        N_population = n_pop,
        Estimate = NA_real_,
        SE = NA_real_,
        Z = NA_real_,
        P = NA_real_,
        IRR = NA_real_,
        Model = "NegativeBinomial",
        Status = "Coefficient unavailable"
      )
    )

  }


  coefficient <- coef_table[
    metric,
    ,
    drop = FALSE
  ]


  ########################################################
  # Return result
  ########################################################

  data.table(
    N = nrow(tmp),
    N_population = n_pop,
    Estimate = coefficient[1, "Estimate"],
    SE = coefficient[1, "Std. Error"],
    Z = coefficient[1, "z value"],
    P = coefficient[1, "Pr(>|z|)"],
    IRR = exp(
      coefficient[1, "Estimate"]
    ),
    Model = "NegativeBinomial",
    Status = "OK"
  )
}


############################################################
# 19. Run population-adjusted models
############################################################

nb_results <- list()

counter <- 1


for (gender_i in levels(sample_df$Gender)) {

  gender_df <- sample_df[
    Gender == gender_i
  ]


  for (metric_i in c(
    "AssemblySize",
    "N50"
  )) {


    for (variant in sv_variants) {

      stat <- run_nb_model(
        gender_df,
        variant,
        metric_i
      )


      stat[
        ,
        `:=`(
          Gender = gender_i,
          Metric = metric_i,
          Variant = variant
        )
      ]


      nb_results[[counter]] <- stat

      counter <- counter + 1
    }
  }
}


nb_dt <- rbindlist(
  nb_results,
  fill = TRUE
)


############################################################
# 20. FDR correction
############################################################
#
# FDR is performed separately within Gender.
#
############################################################

nb_dt[
  ,
  FDR := p.adjust(
    P,
    method = "BH"
  ),
  by = Gender
]


############################################################
# 21. Significance labels
############################################################

nb_dt[
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


############################################################
# 22. Save negative binomial results
############################################################

fwrite(
  nb_dt,
  file.path(
    output_dir,
    "02_population_adjusted_negative_binomial.tsv"
  ),
  sep = "\t"
)


############################################################
# 23. Print negative binomial results
############################################################

cat("\n========================================\n")
cat("Population-adjusted negative binomial models\n")
cat("========================================\n\n")

print(
  nb_dt
)


############################################################
# 24. Extract only INS / DEL results
############################################################

main_nb_results <- nb_dt[
  Variant %in% c(
    "INS",
    "DEL"
  )
]


fwrite(
  main_nb_results,
  file.path(
    output_dir,
    "03_main_INS_DEL_population_adjusted.tsv"
  ),
  sep = "\t"
)


############################################################
# 25. Create human-readable interpretation table
############################################################

interpretation_dt <- copy(
  nb_dt
)


interpretation_dt[
  ,
  Interpretation := fifelse(

    is.na(P),

    "Model unavailable",

    fifelse(

      FDR < 0.05,

      "Assembly quality remains associated with variant count after population adjustment",

      "No significant association after population adjustment"

    )

  )
]


fwrite(
  interpretation_dt,
  file.path(
    output_dir,
    "04_interpretation_table.tsv"
  ),
  sep = "\t"
)


############################################################
# 26. Visualization:
# raw association by Gender and Population
############################################################

plot_long <- melt(
  sample_df,
  id.vars = c(
    "Sample",
    "Gender",
    "Population",
    "N50",
    "AssemblySize"
  ),
  measure.vars = sv_variants,
  variable.name = "Variant",
  value.name = "VariantCount"
)


############################################################
# 27. N50 vs INS / DEL
############################################################

plot_n50_sv <- ggplot(
  plot_long[
    Variant %in% c(
      "INS",
      "DEL"
    )
  ],
  aes(
    x = N50,
    y = VariantCount
  )
) +

  geom_point(
    size = 1.8,
    alpha = 0.65
  ) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.7
  ) +

  facet_grid(
    Variant ~ Gender,
    scales = "free_y"
  ) +

  scale_x_continuous(
    labels = label_number(
      scale = 1e-6,
      suffix = " Mb"
    )
  ) +

  labs(
    x = "Mean haplotype N50",
    y = "SV count"
  ) +

  theme_classic(
    base_size = 15
  ) +

  theme(
    strip.text = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    output_dir,
    "N50_INS_DEL_by_gender.pdf"
  ),
  plot_n50_sv,
  width = 9,
  height = 7
)


############################################################
# 28. AssemblySize vs INS / DEL
############################################################

plot_size_sv <- ggplot(
  plot_long[
    Variant %in% c(
      "INS",
      "DEL"
    )
  ],
  aes(
    x = AssemblySize,
    y = VariantCount
  )
) +

  geom_point(
    size = 1.8,
    alpha = 0.65
  ) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.7
  ) +

  facet_grid(
    Variant ~ Gender,
    scales = "free_y"
  ) +

  scale_x_continuous(
    labels = label_number(
      scale = 1e-9,
      suffix = " Gb"
    )
  ) +

  labs(
    x = "Combined assembly size",
    y = "SV count"
  ) +

  theme_classic(
    base_size = 15
  ) +

  theme(
    strip.text = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold"
    )
  )


ggsave(
  file.path(
    output_dir,
    "AssemblySize_INS_DEL_by_gender.pdf"
  ),
  plot_size_sv,
  width = 9,
  height = 7
)


############################################################
# 29. Model summary for INS and DEL
############################################################

cat("\n========================================\n")
cat("MAIN RESULTS: INS / DEL\n")
cat("========================================\n\n")


print(
  nb_dt[
    Variant %in% c(
      "INS",
      "DEL"
    ),
    .(
      Gender,
      Metric,
      Variant,
      N,
      N_population,
      Estimate,
      IRR,
      P,
      FDR,
      significance,
      Interpretation
    )
  ]
)


############################################################
# 30. End
############################################################

cat("\n========================================\n")
cat("Analysis completed.\n")
cat("========================================\n")

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)

cat(
  "\nMain files:\n",
  "01_descriptive_spearman.tsv\n",
  "02_population_adjusted_negative_binomial.tsv\n",
  "03_main_INS_DEL_population_adjusted.tsv\n",
  "04_interpretation_table.tsv\n"
)

############################################################
# END
############################################################
