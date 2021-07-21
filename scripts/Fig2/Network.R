
path_result <- "results/NATMI_each_patient/Dataframe_AtoE/Dataframe_E.csv"
path_LR <- "results/Fig1/LR_adjPval_meanHR_screened.csv"
path_LR_all <- "results/Fig1/LR_HR_adjPval.csv"
path_outdir <- "results/Fig2/Network"

# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Load library
library(tidyverse)
library(qgraph)
library(ggrepel)

# Count the number of LR pairs
dflr <- read_csv(path_LR)
dflr %>%
  filter(meanHR < 1) %>%
  nrow() -> n_HRL
dflr %>%
  filter(meanHR > 1) %>%
  nrow() -> n_HRH

# Count the number of all LR pairs
dflrall <- read_csv(path_LR_all)
dflrall %>%
  nrow() -> n_LR_all

# Read data
df1 <- read_csv(path_result, col_names = TRUE)

# Calculate Adjusted enrichment
df1 %>%
  mutate(
    adjusted_mean_NormHRL = mean_NormHRL / n_HRL * n_LR_all,
    adjusted_mean_NormHRH = mean_NormHRH / n_HRH * n_LR_all
  ) -> df1

# Save data
write_csv(df1, file.path(path_outdir, "adjusted_mean_enrichment_LRpairs.csv"))


# Network
## Define edges
df1 %>%
  separate(cell_type_pair, sep = "->", into = c("from", "to")) %>%
  mutate(thickness = adjusted_mean_NormHRL) %>%
  select(from, to, thickness) -> Edges_HRL

df1 %>%
  separate(cell_type_pair, sep = "->", into = c("from", "to")) %>%
  mutate(thickness = adjusted_mean_NormHRH) %>%
  select(from, to, thickness) -> Edges_HRH

df1 %>%
  separate(cell_type_pair, sep = "->", into = c("from", "to")) %>%
  mutate(thickness = adjusted_mean_NormHRL - adjusted_mean_NormHRH) %>%
  select(from, to, thickness) -> Edges_diff

# network_cell_type_HRL
pdfname <- file.path(path_outdir, "network_cell_type_HRL.pdf")
pdf(pdfname)
qgraph(Edges_HRL,
  esize = 5,
  theme = "gray", layout = "circle", maximum = max_value,
  title = "LR pairs with HR<1"
)
dev.off()

# network_cell_type_HRH
pdfname <- file.path(path_outdir, "network_cell_type_HRH.pdf")
pdf(pdfname)
qgraph(Edges_HRH,
  esize = 5,
  theme = "gray", layout = "circle", maximum = max_value,
  title = "LR pairs with HR>1"
)
dev.off()

# network_cell_type_HRL
pdfname <- file.path(path_outdir, "network_cell_type_diff.pdf")
pdf(pdfname)
qgraph(Edges_diff,
  esize = 5,
  theme = "gray", layout = "circle",
  title = "Difference of HR<1 from HR>1"
)
dev.off()

# sessionInfo()
sessionInfo()
