############################################################
# APC-P1-specific variants:
# Population distribution / UpSet plot
#
# Reviewer question:
# "Are these variants enriched in specific populations
# that remain underrepresented in existing global references?"
#
# Definition:
# A variant is considered present in a Location if at least
# one individual from that Location carries an ALT allele.
#
# Population is defined by metadata column: Location
############################################################


############################################################
# 0. Libraries
############################################################

library(data.table)
library(ggplot2)
library(ComplexUpset)
library(grid)
library(scales)


############################################################
# 1. Parameters
############################################################

#-----------------------------------------------------------
# Input VCFs
#-----------------------------------------------------------

apc_unique_vcf <-
    "APC_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.vcf.gz"

apc_all_vcf <-
    "APC.chm13v2.norm2.allchr.rmGRCh38.AF.vcf.gz"


#-----------------------------------------------------------
# Metadata
#
# Format:
#
# IID Gender Population Location Region
# CPC0100130 Male Deang Austro-Asiatic China
# ...
#
# IMPORTANT:
# population = Location
#-----------------------------------------------------------

metadata_file <-
    "/path/to/130APC.v3.info.txt"


#-----------------------------------------------------------
# Output
#-----------------------------------------------------------

output_prefix <-
    "APC_P1_specific_population"

output_dir <-
    "APC_population_UPSET"

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


#-----------------------------------------------------------
# bcftools
#-----------------------------------------------------------

bcftools <-
    "/path/to/bcftools"


#-----------------------------------------------------------
# Number of intersections shown
#
# IMPORTANT:
# This does NOT remove populations.
# It only controls how many intersection columns
# are displayed in the UpSet panel.
#-----------------------------------------------------------

n_intersections <- 30


############################################################
# 2. Check input files
############################################################

stop_if_missing <- function(x) {

    if (!file.exists(x)) {

        stop(
            paste0(
                "\nFile does not exist:\n",
                x,
                "\n"
            )
        )
    }

}


stop_if_missing(apc_unique_vcf)
stop_if_missing(apc_all_vcf)
stop_if_missing(metadata_file)


############################################################
# 3. Read sample metadata
############################################################

meta <- fread(
    metadata_file,
    header = TRUE,
    sep = "\t",
    data.table = TRUE
)


############################################################
# 4. Check metadata columns
############################################################

required_meta_columns <- c(
    "IID",
    "Location"
)


missing_meta_columns <-
    setdiff(
        required_meta_columns,
        colnames(meta)
    )


if (length(missing_meta_columns) > 0) {

    stop(
        paste0(
            "Missing metadata columns: ",
            paste(
                missing_meta_columns,
                collapse = ", "
            )
        )
    )

}


############################################################
# 5. Remove duplicated IID
############################################################

if (anyDuplicated(meta$IID)) {

    duplicated_iid <-
        meta[
            duplicated(IID),
            IID
        ]

    warning(
        paste0(
            "Duplicated IID detected:\n",
            paste(
                unique(duplicated_iid),
                collapse = ", "
            )
        )
    )

}


############################################################
# 6. Basic metadata QC
############################################################

if (any(is.na(meta$Location) | meta$Location == "")) {
    stop("Some samples have missing Location.")
}


############################################################
# 7. Population sample size
############################################################

population_size <- meta[
    ,
    .(
        N = uniqueN(IID)
    ),
    by = Location
]


population_size[
    order(-N)
]


############################################################
# 8. Get sample names from APC VCF
############################################################

get_vcf_samples <- function(vcf) {

    cmd <- paste(
        shQuote(bcftools),
        "query -l",
        shQuote(vcf)
    )

    samples <-
        system(
            cmd,
            intern = TRUE
        )

    samples[
        samples != ""
    ]

}


apc_samples <-
    get_vcf_samples(apc_unique_vcf)


############################################################
# 9. Match metadata with VCF samples
############################################################

meta_vcf <- meta[
    IID %in% apc_samples
]


cat(
    "\nNumber of samples in metadata:",
    nrow(meta),
    "\n"
)

cat(
    "Number of samples in APC VCF:",
    length(apc_samples),
    "\n"
)

cat(
    "Number of matched samples:",
    nrow(meta_vcf),
    "\n"
)


############################################################
# 10. Check missing samples
############################################################

missing_from_vcf <-
    setdiff(
        meta$IID,
        apc_samples
    )


if (length(missing_from_vcf) > 0) {

    cat(
        "\nSamples in metadata but not in APC VCF:\n"
    )

    print(
        missing_from_vcf
    )

}


############################################################
# 11. Create population-specific sample lists
############################################################

population_list <-
    split(
        meta_vcf$IID,
        meta_vcf$Location
    )


############################################################
# 12. Population summary
############################################################

population_summary <- meta_vcf[
    ,
    .(
        N = uniqueN(IID)
    ),
    by = Location
]


############################################################
# 13. Extract APC-P1-specific variants
#
# We do NOT read the entire genotype matrix into R.
#
# Instead, for each Location:
#
# 1. subset the VCF to samples from that Location
# 2. retain variants where at least one sample has ALT
# 3. extract variant ID
#
# This is much more memory efficient.
############################################################


cat(
    "\n========================================\n"
)

cat(
    "Extracting population-specific variants\n"
)

cat(
    "========================================\n\n"
)


############################################################
# 14. Create temporary directory
############################################################

tmp_dir <-
    file.path(
        output_dir,
        "tmp_population_vcf"
    )

dir.create(
    tmp_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


############################################################
# 15. Variant ID definition
############################################################
#
# CHROM:POS:REF:ALT
#
# Because the VCF has already been normalized and realigned,
# this is used as the exact variant identity.
#
############################################################

variant_id_format <-
    "%CHROM:%POS:%REF:%ALT"


############################################################
# 16. Loop through populations
############################################################

population_variant_files <-
    character()


for (pop in names(population_list)) {

    cat(
        "Processing:",
        pop,
        "\n"
    )


    ########################################################
    # Population sample list
    ########################################################

    sample_file <-
        file.path(
            tmp_dir,
            paste0(
                "samples_",
                make.names(pop),
                ".txt"
            )
        )


    fwrite(
        data.table(
            IID = population_list[[pop]]
        ),
        sample_file,
        col.names = FALSE
    )


    ########################################################
    # Output variant list
    ########################################################

    variant_file <-
        file.path(
            tmp_dir,
            paste0(
                "variants_",
                make.names(pop),
                ".txt"
            )
        )


    ########################################################
    # bcftools:
    #
    # -S sample list
    #
    # -i 'GT="alt"'
    #
    # This retains a variant if at least one individual
    # in this Location carries an ALT allele.
    #
    ########################################################

    cmd <- paste(
        shQuote(bcftools),
        "query",
        "-S",
        shQuote(sample_file),
        "-i",
        shQuote('GT="alt"'),
        "-f",
        shQuote(
            paste0(
                variant_id_format,
                "\n"
            )
        ),
        shQuote(apc_unique_vcf),
        ">",
        shQuote(variant_file)
    )


    status <-
        system(
            cmd
        )


    if (status != 0) {

        stop(
            paste0(
                "\nbcftools failed for population: ",
                pop
            )
        )

    }


    ########################################################
    # Check number of variants
    ########################################################

    n_variants <-
        length(
            readLines(
                variant_file,
                warn = FALSE
            )
        )


    cat(
        "  samples =",
        length(population_list[[pop]]),
        "\n"
    )

    cat(
        "  variants =",
        n_variants,
        "\n\n"
    )


    population_variant_files[pop] <-
        variant_file

}


############################################################
# 17. Build variant × population membership matrix
############################################################


cat(
    "\nBuilding membership matrix...\n"
)


all_variant_ids <-
    unique(
        unlist(
            lapply(
                population_variant_files,
                function(x) {
                    readLines(
                        x,
                        warn = FALSE
                    )
                }
            )
        )
    )


cat(
    "Total APC-P1-specific variants:",
    length(all_variant_ids),
    "\n"
)

# Total APC-P1-specific variants: 6023905

############################################################
# 18. Initialize matrix
############################################################

membership <-
    matrix(
        FALSE,
        nrow = length(all_variant_ids),
        ncol = length(population_variant_files)
    )


rownames(membership) <-
    all_variant_ids


colnames(membership) <-
    names(population_variant_files)


############################################################
# 19. Fill membership matrix
############################################################

for (pop in names(population_variant_files)) {

    ids <-
        readLines(
            population_variant_files[[pop]],
            warn = FALSE
        )

    ids <-
        unique(
            ids[
                ids != ""
            ]
        )


    membership[
        match(
            ids,
            all_variant_ids
        ),
        pop
    ] <- TRUE

}


############################################################
# 20. Convert to data.table
############################################################

membership_dt <-
    as.data.table(
        membership,
        keep.rownames = "VARIANT"
    )


############################################################
# 21. Remove variants with no population membership
############################################################
#
# This should theoretically be zero because every APC-specific
# variant has at least one ALT allele somewhere in APC.
#
############################################################

membership_dt[
    ,
    total_populations :=
        rowSums(
            .SD
        ),
    .SDcols =
        names(population_variant_files)
]


membership_dt <-
    membership_dt[
        total_populations > 0
    ]


membership_dt[
    ,
    total_populations := NULL
]


############################################################
# 22. Add population size
############################################################

population_summary[
    ,
    N := as.integer(N)
]


############################################################
# 23. Calculate number of variants per population
############################################################

variant_count <- data.table(
    Location =
        names(population_variant_files)
)


variant_count[
    ,
    APC_specific_variants :=
        sapply(
            Location,
            function(x) {
                sum(
                    membership_dt[[x]]
                )
            }
        )
]


############################################################
# 24. Merge population size
############################################################

population_summary <-
    merge(
        population_summary,
        variant_count,
        by = "Location",
        all = TRUE
    )


############################################################
# 25. Fraction of APC-specific variants
############################################################

population_summary[
    ,
    Fraction_APC_specific :=
        APC_specific_variants /
        nrow(membership_dt)
]


############################################################
# 26. Sort populations
############################################################
#
# For visualization:
#
# primary ordering = number of APC-specific variants
# secondary information = population sample size
#
############################################################

population_summary[
    order(
        -APC_specific_variants,
        -N
    )
]


############################################################
# 27. Save population summary
############################################################

fwrite(
    population_summary,
    file.path(
        output_dir,
        paste0(
            output_prefix,
            "_5population_summary.tsv"
        )
    ),
    sep = "\t"
)


############################################################
# 28. Save membership matrix
############################################################

fwrite(
    membership_dt,
    file.path(
        output_dir,
        paste0(
            output_prefix,
            "_5membership_matrix.tsv.gz"
        )
    ),
    sep = "\t"
)

############################################################
# 29. Pre-aggregate data (Vectorized calculation, Avoid loading into memory Bug)
############################################################
library(patchwork)

cat("\nAggregating massive variant data...\n")

intersection_sets <- names(population_variant_files)

membership_dt[
    ,
    (intersection_sets) := lapply(.SD, as.logical),
    .SDcols = intersection_sets
]

# Compress data: Count variants in each intersection
agg_dt <- membership_dt[, .(Count = .N), by = intersection_sets]
agg_dt[, is_empty := rowSums(.SD) == 0, .SDcols = intersection_sets]
agg_dt <- agg_dt[is_empty == FALSE][, is_empty := NULL]


############################################################
# 30. Select the top 30 largest intersections and format
############################################################
agg_dt <- agg_dt[order(-Count)]
n_intersections <- 30
top_agg <- head(agg_dt, n_intersections)

# Assign intersection identifiers ID (Sort in descending order)
top_agg[, Intersect_ID := factor(1:.N, levels = 1:.N)]


############################################################
# 31. Count total variants per population and set order (For the left bar chart)
############################################################
pop_totals <- sapply(intersection_sets, function(pop) {
    sum(agg_dt[get(pop) == TRUE, Count])
})

pop_dt <- data.table(Location = names(pop_totals), Total = pop_totals)
# Sort in ascending order, Because ggplot After flipping axes bars are drawn from bottom to top, Place the largest value at the top
pop_dt <- pop_dt[order(Total)]
pop_levels <- pop_dt$Location


############################################################
# 32. Build the lower membership matrix (Long format)
############################################################
matrix_dt <- melt(
    top_agg,
    id.vars = c("Intersect_ID", "Count"),
    measure.vars = intersection_sets,
    variable.name = "Location",
    value.name = "Present"
)
matrix_dt <- matrix_dt[Present == TRUE]
matrix_dt[, Location := factor(Location, levels = pop_levels)]


############################################################
# 33. Draw the three panels separately (Custom UpSet)
############################################################
cat("Generating Custom UpSet plot via patchwork...\n")

# A. Top bar chart (Top Bar Plot)
p_top <- ggplot(top_agg, aes(x = Intersect_ID, y = Count)) +
    geom_col(fill = "grey35", width = 0.65) +
    geom_text(aes(label = format(Count, big.mark = ",")), vjust = -0.5, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = "Population distribution of APC-P1-specific variants",
        subtitle = paste0("Total distinct APC-P1-specific variants: ", format(nrow(membership_dt), big.mark = ",")),
        y = "Intersection Size",
        x = NULL
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11)
    )

# B. Left bar chart (Left Bar Plot - Use one color per population)
p_left <- ggplot(pop_dt, aes(x = factor(Location, levels = pop_levels), y = Total, fill = Location)) +
    geom_col(width = 0.65) +
    coord_flip() +
    # Flip Y Axis, Draw bars from right to left, Align with the central matrix
    scale_y_reverse(
        expand = expansion(mult = c(0.15, 0)),
        labels = scales::label_comma()
    ) +
    # Place population labels on the right (Place next to the matrix)
    scale_x_discrete(position = "top") +
    labs(y = "Set Size", x = NULL) +
    theme_minimal(base_family = "Arial") +
    theme(
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 11, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 12)
    )

# C. Central membership matrix (Matrix Plot - Use distinct point colors)
# Generate a gray background grid
bg_dt <- as.data.table(expand.grid(
    Intersect_ID = factor(1:n_intersections, levels = 1:n_intersections),
    Location = factor(pop_levels, levels = pop_levels)
))

p_matrix <- ggplot(matrix_dt, aes(x = Intersect_ID, y = Location)) +
    geom_point(data = bg_dt, color = "grey92", size = 3) + #  Gray background points
    geom_line(aes(group = Intersect_ID), color = "grey70", linewidth = 0.6) + #  Connecting lines
    geom_point(aes(color = Location), size = 3.2) + #  Highlighted colored points
    theme_minimal(base_family = "Arial") +
    theme(
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0)
    )


############################################################
# 34. Use cowplot Align and combine panels (Avoid patchwork Bug)
############################################################
# Load cowplot Package
library(cowplot)

# Use plot_grid Arrange four positions (Empty upper-left panel, Upper-right panel, Lower-left panel, Lower-right matrix) Align panels before combining
final_plot <- plot_grid(
    NULL,     p_top,
    p_left,   p_matrix,
    ncol = 2,
    rel_widths = c(0.22, 0.78),  #  Relative panel widths
    rel_heights = c(0.65, 0.35), #  Relative panel heights
    align = "hv",                # h= Align horizontally , v= Align vertically
    axis = "tblr"                #  Align to the top (t),  Bottom (b),  Left (l),  Right (r) Edges
)

############################################################
# 35. Save the final results
############################################################
cat("Saving Custom UpSet Plot...\n")

ggsave(
    filename = file.path(output_dir, paste0(output_prefix, "_5Custom_UpSet.pdf")),
    plot = final_plot,
    width = 16,
    height = 10,
    units = "in",
    device = cairo_pdf
)

ggsave(
    filename = file.path(output_dir, paste0(output_prefix, "_5Custom_UpSet.png")),
    plot = final_plot,
    width = 16,
    height = 10,
    units = "in",
    dpi = 400
)

cat("\n========================================\n")
cat("Analysis fully completed successfully!\n")
cat("========================================\n")
