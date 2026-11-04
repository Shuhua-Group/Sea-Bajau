#!/usr/bin/env Rscript
# Summarize cross-validation errors across ADMIXTURE replicates.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript summarize_admixture_cv.R <outroot>")
outroot <- args[1]
logs <- list.files(file.path(outroot, "02.admixture_runs"), pattern = "\\.log$", recursive = TRUE, full.names = TRUE)
rows <- lapply(logs, function(f) {
  x <- readLines(f, warn = FALSE)
  z <- grep("CV error", x, value = TRUE)
  if (!length(z)) return(NULL)
  m <- regexec("CV error \\(K=([0-9]+)\\):[[:space:]]*([0-9.eE+-]+)", z[length(z)])
  got <- regmatches(z[length(z)], m)[[1]]
  if (length(got) != 3L) return(NULL)
  path <- strsplit(normalizePath(f), "/")[[1]]
  data.frame(Panel = path[length(path) - 3L], K = as.integer(got[2]),
             Replicate = as.integer(sub("rep", "", path[length(path) - 1L])),
             CV_error = as.numeric(got[3]), Log = f)
})
cv <- do.call(rbind, Filter(Negate(is.null), rows))
if (is.null(cv) || !nrow(cv)) stop("No parsable CV errors found.")
cv <- cv[order(cv$Panel, cv$K, cv$Replicate), ]
summary <- do.call(rbind, lapply(split(cv, list(cv$Panel, cv$K), drop = TRUE), function(x) {
  data.frame(Panel = x$Panel[1], K = x$K[1], N_runs = nrow(x),
             CV_median = median(x$CV_error), CV_mean = mean(x$CV_error),
             CV_sd = sd(x$CV_error), CV_min = min(x$CV_error), CV_max = max(x$CV_error),
             CV_q025 = quantile(x$CV_error, .025), CV_q975 = quantile(x$CV_error, .975))
}))
summary <- summary[order(summary$Panel, summary$K), ]
outdir <- file.path(outroot, "03.summary")
write.table(cv, file.path(outdir, "admixture_cv_by_replicate.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(summary, file.path(outdir, "admixture_cv_summary.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
