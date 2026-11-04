# Load required libraries
library(optparse)
library(karyoploteR)

# Parse command line arguments
option_list <- list(
  make_option(c("--wd"), type = "character", help = "Working directory"),
  make_option(c("--prefix"), type = "character", help = "Prefix value"),
  make_option(c("--sample"), type = "character", help = "Sample name"),
  make_option(c("--cyto"), type = "character", help = "Cytobands file path"),
  make_option(c("--gs"), type = "character", help = "Genome structure file path"),
  make_option(c("--readsbam"), type = "character", help = "Reads BAM file path")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Assign parsed arguments to variables
wd <- opt$wd
prefix <- opt$prefix
sample <- opt$sample
cyto <- opt$cyto
gs <- opt$gs
readsbam <- opt$readsbam


# Define a custom BAM density function
kpPlotBAMDensity <- function (karyoplot, data = NULL, window.size = 1e+06, normalize = FALSE,
    ymin = NULL, ymax = NULL, data.panel = 1, r0 = NULL, r1 = NULL,
    col = "gray80", border = NA, clipping = TRUE, ...)
{
    if (!methods::is(karyoplot, "KaryoPlot"))
        stop(paste0("In kpPlotBAMDensity: 'karyoplot' must be a valid 'KaryoPlot' object"))
    if (is.character("data")) {
        data <- Rsamtools::BamFile(file = data)
    }
    if (!methods::is(data, "BamFile"))
        stop(paste0("In kpPlotBAMDensity: 'data' must be a character or a 'BamFile' object."))
    karyoplot$beginKpPlot()
    on.exit(karyoplot$endKpPlot())
    plot.region.lengths <- setNames(width(karyoplot$plot.region),
        as.character(seqnames(karyoplot$plot.region)))
    windows <- tileGenome(plot.region.lengths, tilewidth = window.size,
        cut.last.tile.in.chrom = TRUE)
    seqinfo(windows) <- GenomeInfoDb::Seqinfo(seqnames = seqlevels(windows))
    windows <- shift(windows, shift = start(karyoplot$plot.region[seqnames(windows)]) -
        1)
    dens <- Rsamtools::countBam(data, param = Rsamtools::ScanBamParam(which = windows))$records
    #MODIFIED CODE
    dens <- log(x=dens, base=10)
    #END MODIFIED CODE
    if (normalize == TRUE) {
        total.reads <- sum(Rsamtools::idxstatsBam(file = data)$mapped)
        dens <- dens/total.reads
    }
    if (is.null(ymax)) {
        ymax <- max(dens)
    }
    if (is.null(border) & !is.null(col) & !is.na(col)) {
        border = darker(col, amount = 100)
    }
    if (is.null(col) & !is.null(border) & !is.na(col)) {
        col = lighter(border)
    }
    if (is.na(col) & is.null(border)) {
        border <- "black"
    }
    karyoplot <- kpBars(karyoplot, data = windows, y0 = 0, y1 = dens,
        ymin = ymin, ymax = ymax, data.panel = data.panel, r0 = r0,
        r1 = r1, border = border, col = col, clipping = clipping,
        ...)
    karyoplot$latest.plot <- list(funct = "kpPlotBAMDensity",
        computed.values = list(density = dens, windows = windows,
            max.density = max(dens)))
    invisible(karyoplot)
}

# Define the color palette
palette <- c(
  "#C8CBE8", "#B8EAA2", "#DFED45", "#D77AAD", "#A0E563", "#AC7781", "#EA3FDD",
  "#83B359", "#B5CCA3", "#DBEC81", "#E1A43B", "#E1E8C7", "#E24989", "#60ECBD",
  "#E8E29E", "#E7B1DC", "#C9EAEB", "#A337E5", "#729ACB", "#EBA5A8", "#E65148",
  "#D8C445", "#71EBE7", "#A2E8C6", "#7B3F7F", "#DD845D", "#AF98DE", "#6BE88B",
  "#E5CCC2", "#93C1BD", "#757F5B", "#96909F", "#677DDA", "#63EB3F", "#5EBDDB",
  "#CC61D8", "#645ADB", "#E1B881", "#E18DE4", "#4DAA7F"
)

# Create the PDF
pdf(file = file.path(wd, "11.asm_update/plot_pdf", prefix, sample,
                     paste0(sample, ".alignto_", prefix, ".pdf")), height = 15, width = 12)

# Load cytobands and genome structure
cytobands <- read.table(cyto, header = TRUE)
cytobands <- toGRanges(cytobands)
chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
gs <- toGRanges(gs)

# Plot karyotype
kp <- plotKaryotype(genome = gs, chromosomes = chromosomes, cytobands = cytobands, plot.type = 2)

# Plot haplotype 1
asm1 <- read.table(file.path(wd, "11.asm_update/v2_corstrand", prefix, sample,
                             "bed", paste0(sample, ".hap1.", prefix, ".identifiers_filled.bed")))
asm1_prefix <- paste0(sample, ".hap1.", prefix)
for (chr in chromosomes) {
    asm_subset <- asm1[asm1[, 1] == chr, 1:3]
    asm_subset = toGRanges(asm_subset)
    for (i in seq_along(asm_subset)) {
        color_index <- (i - 1) %% length(palette) + 1
        kpPlotRegions(kp, asm_subset[i], col = palette[color_index], r0 = 0.1, r1 = 0.4, data.panel = 1)
    }
}

asm1_overlap <- read.table(file.path(wd, "11.asm_update/v2_corstrand", prefix, sample,
                                     "bed", paste0(sample, ".hap1.", prefix, ".multi_identifier.bed")))
asm1_overlap = asm1_overlap[,1:3]
asm1_overlap = toGRanges(asm1_overlap)
for (i in seq_along(asm1_overlap)) {
    kpPlotRegions(kp, asm1_overlap[i], col = "red", r0 = 0.1, r1 = 0.4, data.panel = 1)
    }
kpAddLabels(kp, asm1_prefix, r0 = 0.1, r1 = 0.4, data.panel = 1, label.margin = 0.02, cex = 0.3)

# Plot haplotype 2
asm2 <- read.table(file.path(wd, "11.asm_update/v2_corstrand", prefix, sample,
                             "bed", paste0(sample, ".hap2.", prefix, ".identifiers_filled.merged.bed")))
asm2_prefix <- paste0(sample, ".hap2.", prefix)
for (chr in chromosomes) {
    asm_subset <- asm2[asm2[, 1] == chr, 1:3]
    asm_subset = toGRanges(asm_subset)
    for (i in seq_along(asm_subset)) {
        color_index <- (i - 1) %% length(palette) + 1
        kpPlotRegions(kp, asm_subset[i], col = palette[color_index], r0 = 0.1, r1 = 0.4, data.panel = 2)
        }
}
asm2_overlap <- read.table(file.path(wd, "11.asm_update/v2_corstrand", prefix, sample,
                                     "bed", paste0(sample, ".hap2.", prefix, ".multi_identifier.bed")))
asm2_overlap = asm2_overlap[,1:3]
asm2_overlap = toGRanges(asm2_overlap)
for (i in seq_along(asm2_overlap)) {
    kpPlotRegions(kp, asm2_overlap[i], col = "red", r0 = 0.1, r1 = 0.4, data.panel = 2)
    }

kpAddLabels(kp, asm2_prefix, r0 = 0.1, r1 = 0.4, data.panel = 2, label.margin = 0.02, cex = 0.3)

# Plot reads density
kp <- kpPlotBAMDensity(kp, data = readsbam, window.size = 1e5, col = "orchid", r0 = 0.6, r1 = 1, data.panel = 1)
kpAxis(kp, ymax = kp$latest.plot$computed.values$max.density, r0 = 0.6, r1 = 1, cex = 0.2)
# Add legend
legend(x = "bottomright", fill = c("orchid", "red"), legend = c("lg(read_density)", "duplicated_identifier"))
# Close PDF
dev.off()
