# Summarize population variation and demographic history.

library(relater)
library(ggplot2)
library(dplyr)
library(patchwork)   #  For  4×5  Order

# ---- Step 1. Read and organize data ----
coal <- read.coal("/path/to/allsample.rm2018CellAFR.N30000.chr1_22.popsize.seed1.pairwise.coal")

coal$popsize <- 0.5 / coal$haploid.coalescence.rate
coal$epoch.start <- 28 * coal$epoch.start

coal <- coal[coal$epoch.start > 0 & coal$haploid.coalescence.rate > 0, ]
coal <- subset(coal, group1 == group2)

# ---- Step 2. Custom colors ----
my_colors <- c(
  "2018Cell.Bajau" = "#4197d8",
  "36Bajau"        = "#4169E1",
  "2018Cell.Saluan"= "#6D3F95",
  "25Malay"        = "#fc8d59",
  "EGDP.Borneo"    = "#984ea3",
  "EGDP.Luzon"     = "#bebada",
  "EGDP.Myanmar"   = "#B3B3B3",
  "EGDP.Vietnam"   = "#737373",
  "ho1926.Mbuti"   = "#B99B60",
  "ho1926.Yoruba"  = "#b15928",
  "Pacific.Barito" = "#d53e4f",
  "Pacific.Besemah"= "#FFB6C1",
  "Pacific.Kalimantan"="#e78ac3",
  "Pacific.Mamanwa"="#9e70a9",
  "Pacific.Manobo" = "#b8a3e1",
  "Pacific.Semende"="#FF69B4",
  "PanMLS.Jakun"   = "#1aac7c",
  "PanMLS.Temiar"  = "#62C15D"
)

# ---- Step 3. Set target populations ----
target <- "36Bajau"

# ---- Step 4. Generate 17 Figures, Each figure = 36Bajau vs Another population ----
all_groups <- sort(unique(coal$group2))
other_groups <- setdiff(all_groups, target)

plots <- list()

for (g in other_groups) {
  df_plot <- coal %>% filter(group2 %in% c(target, g))

  p <- ggplot(df_plot) +
    geom_step(aes(x = epoch.start, y = popsize, color = group2), linewidth = 0.4) +
    geom_vline(xintercept = c(2000, 3000, 4000,5000),
               linetype = "dashed",
               linewidth = 0.3,
               color = "grey40") +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = my_colors) +
    labs(title = paste0(target, " vs ", g),
         x = "Year",
         y = "Effective Population Size (Ne)") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold"))

  plots[[g]] <- p
}

# ---- Step 5. Use patchwork Combine into 4×5(20 Grid) ----
final_plot <- wrap_plots(plots, ncol = 4, nrow = 5)

# Show
final_plot

# Save
ggsave("RELATE.Ne.36Bajau_all.retry.pdf", final_plot, width = 14, height = 18)
