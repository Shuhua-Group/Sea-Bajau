# Process variant calls and summarize variant catalogs.

library(data.table)
library(ggplot2)
library(dplyr)
# library(ggdist) # For plotting half-violin (Raincloud Core)
# library(patchwork)
library(cowplot)

# 1. Read data
df <- fread("APC_uniq.norm2.chm13v2.allchr.individual_burden.tsv")

# 2. Main filtering step: Exclude N < 10 populations
pop_counts <- df %>%
    group_by(Location) %>%
    summarise(N = n_distinct(IID))

keep_pops <- pop_counts %>%
    filter(N >= 10) %>%
    pull(Location)

df_filtered <- df %>%
    filter(Location %in% keep_pops)

cat(sprintf("Retained %d populations with N >= 10\n", length(keep_pops)))

# 3. Choose the variables to plot Class Order
class_order <- c("COMPLEX", "DEL", "INS", "INDEL", "MNP", "SNP")

# 4. Define one Class plotting function (Pure ggplot2 Native version)
plot_raincloud <- function(cls_name) {

    # Extract the current Class data
    tmp_df <- df_filtered %>% filter(CLASS == cls_name)

    # Based on the current Class Of Burden Order populations by decreasing median
    median_order <- tmp_df %>%
        group_by(Location) %>%
        summarise(Median_Burden = median(Burden)) %>%
        arrange(Median_Burden) %>%
        pull(Location)

    tmp_df$Location <- factor(tmp_df$Location, levels = median_order)

    # Use native ggplot2 Plot
    p <- ggplot(tmp_df, aes(x = Burden, y = Location, fill = Location)) +

        # A. Kernel density at the bottom (Full violin plot)
        geom_violin(
            width = 0.8,
            alpha = 0.4,
            color = NA,
            trim = FALSE
        ) +

        # B. Individual observations in the center (Jitter)
        geom_jitter(
            aes(color = Location),
            height = 0.15,
            size = 1.2,
            alpha = 0.6
        ) +

        # C. Box plot at the top (Boxplot)
        geom_boxplot(
            width = 0.15,
            fill = "white",
            alpha = 0.8,
            outlier.shape = NA,
            color = "black"
        ) +

        labs(
            title = cls_name,
            x = if(cls_name %in% c("MNP", "SNP")) "Number of variants per individual" else NULL,
            y = NULL
        ) +
        scale_x_continuous(labels = scales::comma) +
        theme_classic(base_size = 14) +
        theme(
            legend.position = "none",
            plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.y = element_text(color = "black", face = "bold")
        )

    return(p)
}

# 5. Generate all 6 Records Class plots
plot_list <- list()
for (cls in class_order) {
    if (cls %in% unique(df_filtered$CLASS)) {
        plot_list[[cls]] <- plot_raincloud(cls)
    }
}

# 6. [ Main update ] Use cowplot Combine panels and save
# plotlist Accept a list, ncol = 2 Two columns
# align = "hv" Force horizontal (h) and vertical (v) Align plot panels
final_plot <- plot_grid(
    plotlist = plot_list,
    ncol = 2,
    align = "hv",
    axis = "tblr"
)

ggsave(
    "APC_uniq.Figure_Panel_C_Individual_Burden.pdf",
    final_plot,
    width = 14,
    height = 12,
    device = cairo_pdf
)
