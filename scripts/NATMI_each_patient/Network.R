
path_result <- "results/NATMI_each_patient/Dataframe_AtoE/Dataframe_E.csv"
path_LR <- "results/Fig1/LR_adjPval_meanHR_screened.csv"
path_LR_all <- "results/Fig1/LR_HR_adjPval.csv"
path_outdir <- "results/NATMI_each_patient/Network"

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

# Scatter plot
min_value <- min(c(df1$adjusted_mean_NormHRL, df1$adjusted_mean_NormHRH))
max_value <- max(c(df1$adjusted_mean_NormHRL, df1$adjusted_mean_NormHRH))

df1 %>%
  ggplot(aes(adjusted_mean_NormHRL, adjusted_mean_NormHRH, label = cell_type_pair)) +
  geom_point(aes(color = adjusted_mean_NormHRL / adjusted_mean_NormHRH > 2)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0, slope = 2, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1 / 2, linetype = "dashed") +
  geom_text_repel(data = df1 %>% filter(adjusted_mean_NormHRH / adjusted_mean_NormHRL > 2 | adjusted_mean_NormHRL / adjusted_mean_NormHRH > 2)) +
  theme(legend.position = "bottom") +
  lims(x = c(min_value, max_value), y = c(min_value, max_value)) +
  labs(
    title = "Adjusted mean enrichment of LR pairs",
    x = "LR pairs with HR<1", y = "LR pairs with HR>1"
  ) -> g1
g1
ggsave(file.path(path_outdir, "scatter_adjusted_mean_enrichment_HRL_vs_HRH.pdf"), g1)


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

df1 %>%
  separate(cell_type_pair, sep = "->", into = c("from", "to")) %>%
  mutate(thickness = log2(adjusted_mean_NormHRL + 1e-5) - log2(adjusted_mean_NormHRH + 1e-5)) %>%
  select(from, to, thickness) -> Edges_fc

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

# network_cell_type_HRL
pdfname <- file.path(path_outdir, "network_cell_type_fc.pdf")
pdf(pdfname)
qgraph(Edges_fc,
  esize = 5,
  theme = "gray", layout = "circle",
  title = "log2 fold change of HR<1 from HR>1"
)
dev.off()

# sessionInfo()
sessionInfo()
