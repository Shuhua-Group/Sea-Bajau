############################################################
# Variant class AFS plot
############################################################


library(data.table)
library(ggplot2)
library(scales)


############################################################
# parameters
############################################################

# unique_file <- "APC_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.filter208.tsv"
# all_file <- "APC.chm13v2.norm2.allchr.rmGRCh38.AF.filter208.tsv"
# output_prefix <- "APC"
# sample_number <- 130

# unique_file <- "CPC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.filter370.tsv"
# all_file <- "CPC2.chm13v2.norm2.allchr.rmGRCh38.AF.filter370.tsv"
# output_prefix <- "CPC2"
# sample_number <- 231

unique_file <- "HPRC2_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.filter370.tsv"
all_file <- "HPRC2.chm13v2.norm2.allchr.rmGRCh38.AF.filter370.tsv"
output_prefix <- "HPRC2"
sample_number <- 231

############################################################
# read data
############################################################

unique <- fread(unique_file)
all <- fread(all_file)

############################################################
# add group
############################################################


unique[, group := "Graph specific"]

all[, group := "All"]


# Ensure numeric columns use numeric
num_cols <- c(
    "AC",
    "AN",
    "AF",
    "REF_LEN",
    "ALT_LEN",
    "LEN_DIFF"
)

unique[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
all[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

############################################################
# merge
############################################################

df <- rbind(
    unique,
    all,
    fill=TRUE
)

############################################################
# filter
############################################################

df <- df[
    AF > 0 &
    AF <=1
]

############################################################
# variant class order
############################################################

class_order <- c(
    "COMPLEX",
    "DEL",
    "INS",
    "INDEL",
    "MNP",
    "SNP"
)


df[, CLASS := factor(
    CLASS,
    levels=class_order
)]


############################################################
# AFS bins
############################################################

af_breaks <- c(
    0,
    0.005,
    0.01,
    0.02,
    0.05,
    0.10,
    0.20,
    0.50,
    1
)

af_labels <- c(
    "0–0.005",
    "0.005–0.01",
    "0.01–0.02",
    "0.02–0.05",
    "0.05–0.10",
    "0.10–0.20",
    "0.20–0.50",
    "0.50–1.00"
)


df[
    ,
    AFBIN := cut(
        AF,
        breaks = af_breaks,
        labels = af_labels,
        include.lowest = TRUE,
        right = TRUE
    )
]


df[, AFBIN := factor(
    AFBIN,
    levels = af_labels
)]


############################################################
# count
############################################################


afs <- df[
    !is.na(AFBIN),
    .(N = .N),
    by = .(CLASS, group, AFBIN)
]

totals <- afs[
    ,
    .(Total = sum(N)),
    by = .(CLASS, group)
]

afs <- merge(
    afs,
    totals,
    by = c("CLASS", "group"),
    all.x = TRUE
)

afs[
    ,
    Fraction := N / Total
]


############################################################
# plot
############################################################


p <- ggplot(
    afs,
    aes(
        x=AFBIN,
        y=Fraction,
        color=group,
        group=group
    )
)+


geom_line(
    linewidth=1
)+


geom_point(
    size=3
)+


scale_y_log10(
    labels = scales::label_scientific()
)+


facet_wrap(
    ~CLASS,
    ncol=2
)+


labs(
    x="Allele frequency spectrum",
    y="Fraction of variants (log10)",
    color=NULL
)+


theme_classic(
    base_size=18
)+


theme(
    strip.text =
        element_text(
            face="bold",
            size=20
        ),

    axis.text.x =
        element_text(
            angle=45,
            hjust=1
        ),

    legend.position="top"
)



############################################################
# save
############################################################


ggsave(
    paste0(output_prefix,
           "_AFS.pdf"),
    p,
    width=12,
    height=10
)
