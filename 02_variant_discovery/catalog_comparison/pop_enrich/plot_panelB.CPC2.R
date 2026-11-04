############################################################
# 1. Load dependencies
############################################################
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

cat("1.  Reading and filtering data ...\n")

############################################################
# 2. Read and filter data
############################################################
df <- fread("CPC2_uniq.allchr.AF.miss0.2.tsv")
variant_class <- fread("../CPC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.tsv")

# Keep reliable estimates only AF Estimate (Check case sensitivity)
df <- df[PASS == "True"]

# 1. Coerce N columns to integers (Integer)
df[, N := as.integer(N)]
df <- df[N >= 10]

# Inspect available entries Population Type
unique(df$Population)


############################################################
# 3. Clean data and validate joins
############################################################
cat("2.  Format and merge variant annotations ...\n")

# Coerce to data.frame Avoid issues in underlying C stack handling
df <- as.data.frame(df)
variant_class <- as.data.frame(variant_class)

# Unify POS as integer
df$POS <- as.integer(df$POS)
variant_class$POS <- as.integer(variant_class$POS)

# Unify CHROM format and trim whitespace (Validate and clean data)
df$CHROM <- str_trim(gsub("^chr", "", df$CHROM, ignore.case = TRUE))
variant_class$CHROM <- str_trim(gsub("^chr", "", variant_class$CHROM, ignore.case = TRUE))
df$REF <- str_trim(df$REF)
df$ALT <- str_trim(df$ALT)
variant_class$REF <- str_trim(variant_class$REF)
variant_class$ALT <- str_trim(variant_class$ALT)

# Efficient validated join
df <- inner_join(
    df,
    variant_class,
    by = c("CHROM", "POS", "REF", "ALT")
)

# After merging convert back to data.table for efficient summaries
setDT(df)

# Remove N<10 rows


############################################################
# 4. Set variant-class order
############################################################
cat("3.  Summarize variant frequencies by class ...\n")

class_order <- c("COMPLEX", "DEL", "INS", "INDEL", "MNP", "SNP")

df[, CLASS := factor(CLASS, levels = class_order)]

# If there are unknown CLASS will become NA, Exclude these entries
df <- df[!is.na(CLASS)]


############################################################
# 5. Calculate counts and proportions
############################################################
# Step 1: Count records for each combination (N)
plot_df <- df[, .(N = .N), by = .(CLASS, Population, AF_category)]

# Step 2: Calculate variant-class proportions within each population (Fraction)
plot_df[, Fraction := N / sum(N), by = .(CLASS, Population)]


############################################################
# 6. Draw a percentage stacked bar chart
############################################################
cat("4.  Generating and saving figures ...\n")
# Define AF_category fixed order
af_order <- c("Absent", "Rare", "Low-frequency", "Private-common")

# Use AF_category Convert to an ordered factor
plot_df[, AF_category := factor(AF_category, levels = af_order)]

p <- ggplot(plot_df, aes(x = Population, y = Fraction, fill = AF_category)) +
    geom_col(position = "fill", width = 0.8) +
    facet_wrap(~CLASS, ncol = 2) +

    # Use Y Format axis labels as percentages
    scale_y_continuous(labels = percent) +

    # Use a publication palette (Set2 Or Set1 Etc.)
    scale_fill_brewer(palette = "Set2") +

    labs(
        x = NULL,
        y = "Percentage of CPC2-specific variants",
        fill = "Population AF"
    ) +

    theme_classic(base_size = 16) +
    theme(
        strip.text = element_text(face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right"
    )

############################################################
# 7. Save results
############################################################
ggsave(
    filename = "CPC2_unique_population_AF_category.pdf",
    plot = p,
    width = 14,
    height = 10,
    device = cairo_pdf  #  Use  cairo_pdf  Preserve vector text and lines in layout software
)

p_line <- ggplot(plot_df, aes(x = Population, y = Fraction, color = AF_category, group = AF_category)) +
    # Draw lines, Increase the line width
    geom_line(linewidth = 1.2, alpha = 0.8) +
    # Draw data points, Increase marker size, Use larger markers
    geom_point(size = 3) +

    facet_wrap(~CLASS, ncol = 2) +
    scale_y_continuous(labels = percent) +
    scale_color_brewer(palette = "Set2") +

    labs(
        x = NULL,
        y = "Percentage of CPC2-specific variants",
        color = "Population AF"
    ) +

    theme_classic(base_size = 16) +
    theme(
        strip.text = element_text(face = "bold", size = 16),
        # Format the line plot, X Population names on the axis: 45 Tilt labels by this angle for readability
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right",
        # Add light horizontal grid lines, Make comparisons easier Y Show percentages on the axis
        panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
    )

ggsave("CPC2_unique_population_AF_category_LINE.pdf", p_line, width = 14, height = 10, device = cairo_pdf)

# Use Absent Remove these entries before plotting
############################################################
# 8. Extended version: Exclude "absent" Redraw bars after exclusion
############################################################
cat("5.  Calculating bars after excluding  absent  bar chart after filtering ...\n")

# 1. Exclude AF_category As "absent" rows
plot_df_no_absent <- plot_df[AF_category != "Absent"]

# 2. Recalculate after excluding absent within-group percentages after filtering
# For each Population within rare + Low-frequency + Private-common will sum to 100%
plot_df_no_absent[, Fraction := N / sum(N), by = .(CLASS, Population)]

# 3. Draw the updated stacked bar chart
p_bar_no_absent <- ggplot(
    plot_df_no_absent,
    aes(x = Population, y = Fraction, fill = AF_category)
) +
    geom_col(position = "fill", width = 0.8) +
    facet_wrap(~CLASS, ncol = 2) +

    scale_y_continuous(labels = percent) +
    scale_fill_brewer(palette = "Set2") +

    labs(
        x = NULL,
        y = "Percentage of Present Variants (excluding absent)", #  Updated  Y  Use precise axis labels
        fill = "Population AF"
    ) +

    theme_classic(base_size = 16) +
    theme(
        strip.text = element_text(face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "right"
    )

############################################################
# 4. Save the updated figure
############################################################
ggsave(
    filename = "CPC2_unique_population_AF_category_NO_ABSENT.pdf",
    plot = p_bar_no_absent,
    width = 14,
    height = 10,
    device = cairo_pdf
)

cat(" Exclude  Absent  figure saved to  CPC2_unique_population_AF_category_NO_ABSENT.pdf\n")
