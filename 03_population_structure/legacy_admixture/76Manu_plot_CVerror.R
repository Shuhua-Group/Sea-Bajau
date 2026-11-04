# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)

# ==========================================
# 1. Read and process data
# ==========================================
# Expected input location: cv_error_results.txt In, with a header (directory, K_value, cv_error)
cv_df <- read_tsv("/path/to/76Bajau_CVerror.txt")


# Calculate for each K mean CV error and standard deviation (SD)
cv_summary <- cv_df %>%
  group_by(K_value) %>%
  summarise(
    mean_cv = mean(cv_error),
    sd_cv = sd(cv_error),
    .groups = 'drop'
  )

# Find the mean automatically CV error Best minimum K Value
optimal_k <- cv_summary$K_value[which.min(cv_summary$mean_cv)]

# ==========================================
# 2. Plot with error bars CV Error Curve
# ==========================================
p_cv_errorbar <- ggplot(cv_summary, aes(x = K_value, y = mean_cv)) +

  # 1. Highlight the best value with a vertical dashed line: K Value (Draw behind other elements, Avoid obscuring data points)
  geom_vline(xintercept = optimal_k, linetype = "dashed",
             color = "#4DBBD5", linewidth = 1) +

  # 2. Add Error bar (Mean ± Standard deviation)
  # Set an appropriate width (width) Use narrow horizontal caps, Use dark gray
  geom_errorbar(aes(ymin = mean_cv - sd_cv, ymax = mean_cv + sd_cv),
                width = 0.3, color = "gray30", linewidth = 0.8) +

  # 3. Connect the means
  geom_line(color = "#E64B35", linewidth = 1.2) +

  # 4. Draw mean values (Draw in the foreground)
  geom_point(color = "#E64B35", size = 3) +

  # Format axes and labels
  scale_x_continuous(breaks = seq(2, 20, by = 2)) +
  labs(
    x = "Number of Ancestral Populations (K)",
    y = "Cross-Validation (CV) Error",
    title = "ADMIXTURE Cross-Validation Error",
    subtitle = paste0("Optimal K = ", optimal_k, " (Mean \u00B1 SD across 10 seeds)")
  ) +

  # Use a classic theme
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray30", face = "italic", size = 12)
  )

# Print and save
print(p_cv_errorbar)
ggsave("Figure_ADMIXTURE_CV_ErrorBar.pdf", plot = p_cv_errorbar, width = 20, height = 4, dpi = 300)

