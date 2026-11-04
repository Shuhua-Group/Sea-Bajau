# First, Identify the optimal m

# Directory path
path_llik <- "/path/to/TreeMix"

files <- list.files(path_llik, pattern = "\\.llik$", full.names = TRUE)

ll_data <- data.frame(file = character(), m = integer(), LL = numeric(), stringsAsFactors = FALSE)

for(f in files){
  lines <- readLines(f)
  # Found Exiting Row
  exit_line <- lines[grepl("^Exiting ln\\(likelihood\\)", lines)]
  if(length(exit_line) == 1){
    # Extract m
    m <- as.integer(gsub("Exiting ln\\(likelihood\\) with (\\d+) migration events:.*", "\\1", exit_line))
    # Extract LL
    LL <- as.numeric(gsub("Exiting ln\\(likelihood\\) with \\d+ migration events:\\s*", "", exit_line))

    ll_data <- rbind(ll_data, data.frame(file = basename(f), m = m, LL = LL))
  }
}

ll_data






setwd("/path/to/TreeMix") # of course this needs to be adjusted
prefix="4panel.chrA.shapeit5"





library(RColorBrewer)
library(R.utils)
source("/path/to/plotting_funcs.R") # here you need to add the path

par(mfrow=c(2,3))
for(edge in 0:5){
  pdf(paste0(prefix,".",edge,".pdf"),width=6,height=6)
  plot_tree(cex=0.8,paste0(prefix,".",edge))
  title(paste(edge,"edges"))
  dev.off()
}

