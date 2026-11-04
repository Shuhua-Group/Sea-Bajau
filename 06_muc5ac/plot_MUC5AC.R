# -------------------------------
# Load packages
# -------------------------------
library(gggenomes)
library(tidyverse)
library(ggforce)

setwd("/path/to/07.statistics")

# -------------------------------
# 1. Read HLA Gene metadata
# genes The table must contain at least: seq_id, start, end, strand, name, type
# -------------------------------
genes <- read_tsv("MUC5AC.CDS.gene.tsv",
                  col_types = cols()) %>%
  mutate(
    start = as.integer(start),
    end   = as.integer(end),
    type  = "CDS"   # gggenomes  Requires  type  Fields
  )

# -------------------------------
# 2. Read HLA Sequence metadata
# seqs The table must contain at least: seq_id, seq_desc, length
# -------------------------------
seqs <- read_tsv("MUC5AC.seqs.tsv",
                 col_types = cols()) %>%
  mutate(
    seq_desc = as.character(seq_desc),
    seq_desc = ifelse(is.na(seq_desc), "NA", seq_desc),
    length = as.integer(length)
  )

# -------------------------------
# 3. Read BLAST Results, Generate a link table
# -------------------------------
blast_file <- "MUC5AC.blast.out"

blast <- read_tsv(
  blast_file,
  col_names = c(
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    qseqid = sub("#.*", "", qseqid),
    sseqid = sub("#.*", "", sseqid),
    direction = ifelse(send > sstart, 1, -1)
  )

links <- blast %>%
  transmute(
    seq_id = qseqid,
    seq_id2 = sseqid,
    start = qstart,
    end = qend,
    start2 = sstart,
    end2 = send,
    pident = pident,
    direction = direction
  )
# # Filter records
# links <- links %>%
# filter(pident >= 98)

# -------------------------------
# 3. Read meth Information
# -------------------------------
SV_file<- read_tsv("MUC5AC.SV.bed.tsv",
                   col_names = c("seq_id", "start", "end", "score"),
                   skip = 0, #  If the original file has an incorrect header ,  keep this line ;  If there is no header ,  Change to  skip = 0  or remove this line
                   show_col_types = FALSE) %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end),
    score = as.numeric(score)
  )

# -------------------------------
# 4. identity Discretize into 4 colors
# -------------------------------
identity_colors <- c("#ABCBD8", "#CDDBDB", "#DFE2C3", "#ECDBA3", "#F0C090", "#F19994")
# identity_colors <- c("#ECDBA3", "#F0C090", "#F19994")



links <- links %>%
  mutate(
    pident_group = as.character(cut(
      pident,
      breaks = c(88, 90, 92,94,96,98, 100),
      include.lowest = TRUE,
      labels = identity_colors
    ))
  )

# -------------------------------
# 5. For genes Set colors (Based on name)
# -------------------------------
gene_colors <- c(
  "MUC5AC" = "#e76f51"
)

# Add fill Column
genes <- genes %>%
  mutate(fill = gene_colors[name])

# -------------------------------
# 6. Plot
# -------------------------------
p1 <- gggenomes(
  genes = genes,
  seqs = seqs,
  links = links,
  feats = list(SV=SV_file)
) +
  geom_seq() +
  geom_feat(aes(z = score, linetype = "SV annotation"),
              geom = "linerange", feats(SV),
              position=position_nudge(y= -0.1),
              height = .05) +
  geom_bin_label() +
  geom_gene(aes(fill = fill),position = position_nudge(y= .0), size = 1) + # gene fill
  geom_gene_tag(aes(label=name), nudge_y=0.1, check_overlap = TRUE) +
  geom_link_curved(aes(fill = pident_group, color = pident_group),       # link color
                   curve = 10, alpha = 0.5) +
  scale_fill_identity() +   # gene  Of  fill  Use directly  hex
  scale_color_identity() + # link  use the existing colors  hex
  theme(
    legend.position = "right"
  )
p1

# -------------------------------
# 7. Overview plot
# -------------------------------
pdf("GSTM.hap_align.pdf", height = 6, width = 20)
p1
dev.off()

