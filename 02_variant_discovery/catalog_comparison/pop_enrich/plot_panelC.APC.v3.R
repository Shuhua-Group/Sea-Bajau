# Process variant calls and summarize variant catalogs.

library(data.table)
library(ggplot2)
library(dplyr)
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

# ==========================================================
# Unify Population sort order
# ==========================================================
pop_order_df <- df_filtered %>%
    group_by(IID, Location) %>%
    summarise(Total_Burden = sum(Burden), .groups = 'drop') %>%
    group_by(Location) %>%
    summarise(Median_Total = median(Total_Burden)) %>%
    arrange(Median_Total)

global_pop_order <- pop_order_df$Location

# Use the same value across panels: Location Convert to a common fixed-order Factor
df_filtered$Location <- factor(df_filtered$Location, levels = global_pop_order)

# 3. Choose the variables to plot Class Order
class_order <- c("COMPLEX", "DEL", "INS", "INDEL", "MNP", "SNP")

# 4. Define one Class plotting function
plot_raincloud <- function(cls_name) {

    # Extract the current Class data
    tmp_df <- df_filtered %>% filter(CLASS == cls_name)

    # ==========================================================
    # [ Add ]: Calculate globally Kruskal-Wallis Test P-value
    # ==========================================================
    kw_test <- kruskal.test(Burden ~ Location, data = tmp_df)
    pval <- kw_test$p.value

    # Format P-value Labels
    if (pval < 0.001) {
        p_label <- "Kruskal-Wallis, p < 0.001"
    } else {
        p_label <- sprintf("Kruskal-Wallis, p = %.3f", pval)
    }

    # Calculate the extreme-value threshold (99.5% Quantile) Truncate X Axis
    upper_limit <- quantile(tmp_df$Burden, 0.995)

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

        coord_cartesian(xlim = c(min(tmp_df$Burden), upper_limit)) +

        labs(
            title = cls_name,
            subtitle = p_label, #  Use  P-value  Add as a subtitle
            x = if(cls_name %in% c("MNP", "SNP")) "Number of variants per individual" else NULL,
            y = NULL
        ) +
        scale_x_continuous(labels = scales::comma) +
        theme_classic(base_size = 14) +
        theme(
            legend.position = "none",
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(face = "italic", size = 11, hjust = 0.5, color = "grey40"), #  Use gray italic text for the subtitle ,  Use consistent formatting
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

# 6. Use cowplot Combine panels and save
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

cat("Success! Generated plots with global Kruskal-Wallis p-values.\n")
