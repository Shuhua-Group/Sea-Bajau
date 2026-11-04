############################################################
# Variant class AFS plot
############################################################


library(data.table)
library(ggplot2)
library(scales)


############################################################
# parameters
############################################################

unique_file <- "APC_uniq.norm2.chm13v2.allchr.rmGRCh38.AF.filter208.tsv"
all_file <- "APC.chm13v2.norm2.allchr.rmGRCh38.AF.filter208.tsv"
output_prefix <- "APC"
sample_number <- 130

############################################################
# read data
############################################################

unique <- fread(unique_file)
all <- fread(all_file)

############################################################
# add group
############################################################


unique[, group := "APC-P1-specific"]

all[, group := "All APC"]



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


df[, AFBIN := NA_character_]

# singleton
df[
    AC == 1,
    AFBIN := "Singleton"
]


# doubleton
df[
    AC == 2,
    AFBIN := "Doubleton"
]


# remaining bins

df[
    is.na(AFBIN) &
    AF >0.01 &
    AF <=0.02,
    AFBIN := "0.01-0.02"
]


df[
    is.na(AFBIN) &
    AF >0.02 &
    AF <=0.05,
    AFBIN := "0.02-0.05"
]


df[
    is.na(AFBIN) &
    AF >0.05 &
    AF <=0.1,
    AFBIN := "0.05-0.1"
]


df[
    is.na(AFBIN) &
    AF >0.1 &
    AF <=0.2,
    AFBIN := "0.1-0.2"
]


df[
    is.na(AFBIN) &
    AF >0.2 &
    AF <=0.5,
    AFBIN := "0.2-0.5"
]


df[
    is.na(AFBIN) &
    AF >0.5,
    AFBIN := "0.5-1"
]



############################################################
# bin order
############################################################


bin_order <- c(
    "Singleton",
    "Doubleton",
    "0.01-0.02",
    "0.02-0.05",
    "0.05-0.1",
    "0.1-0.2",
    "0.2-0.5",
    "0.5-1"
)


df[, AFBIN := factor(
    AFBIN,
    levels=bin_order
)]



############################################################
# count
############################################################


afs <- df[
    !is.na(AFBIN),
    .(N=.N
    ),
    by=.(CLASS,group,AFBIN)
]


totals <- afs[
    ,
    .(Total = sum(N)),
    by = .(CLASS, group)
]

afs <- merge(
    afs,
    totals,
    by=c("CLASS","group")
)


afs[
    ,
    Fraction:=N/Total
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
    labels=scientific
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
