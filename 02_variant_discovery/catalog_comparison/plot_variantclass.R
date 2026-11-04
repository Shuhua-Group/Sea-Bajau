# Process variant calls and summarize variant catalogs.

library(ggplot2)
library(dplyr)



# df <- read.table(
# "APC_variant_class.rmGRCh38.tsv.fraction.tsv",
# header=T,
# sep="\t"
#)
# df <- read.table(
# "CPC2_variant_class.rmGRCh38.tsv.fraction.tsv",
# header=T,
# sep="\t"
#)
df <- read.table(
"HPRC2_variant_class.rmGRCh38.tsv.fraction.tsv",
header=T,
sep="\t"
)


# Set variant class Order
class_order <- c(
    "COMPLEX",
    "DEL",
    "INS",
    "INDEL",
    "MNP",
    "SNP"
)


df$CLASS <- factor(
    df$CLASS,
    levels = class_order
)


pdf("HPRC2_variant_class.rmGRCh38.pdf")

ggplot(
    df,
    aes(
        x = DATASET,
        y = FRACTION,
        fill = CLASS
    )
) +
    geom_bar(
        stat = "identity",
        position = "stack"
    ) +
    theme_classic()
dev.off()
