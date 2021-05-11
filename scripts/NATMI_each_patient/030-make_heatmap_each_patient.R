# This script takes `ExtractEdges.py` results as input and generates heatmap PDF as output.

# Definition of input and output paths
path_natmi_result_base <- "results/results/NATMI_each_patient/ExtractEdges"
path_outdir <- "results/NATMI_each_patient/heatmap"

# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Load library
library(tidyverse)

path_natmi_results <- list.files(path_natmi_result_base, full.names = TRUE)

# for each patient
for (path_natmi_result in path_natmi_results) {

  # List up files
  inputfiles <- list.files(path_natmi_results)

  # Read data
  # TODO
  lapply(inputfiles, function(x) {
    read_tsv(x, col_names = TRUE) %>%
      mutate(XX = YYY)
  }) %>% bind_rows() -> df1

  # Make plot
  # TODO
  df1 %>%
    ggplot(aes(hoge, fuga, z = XXXX)) +
    geom_tile() +
    labs(title = sampleid) -> g

  # Save plot
  path_pdf <- file.path(path_outdir, paste0(sampleid, "_ExtractEdges_heatmap.pdf"))
  ggsave(g, path_pdf)

  # Save tsv for the plot
  path_tsv <- file.path(path_outdir, paste0(sampleid, "_ExtractEdges_heatmap.tsv"))
  df1 %>%
    write_tsv(path_tsv)
}

# sessioninfo
sessionInfo()
