
path_result <- "results/NATMI_each_patient/summary_number_of_patient.csv"
path_outdir <- "results/NATMI_each_patient/Network_singleLR"

# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Load library
library(tidyverse)
library(ggrepel)
library(qgraph)

# Read data
df1 <- read_csv(path_result) %>%
  replace_na(list(HRL = 0, HRH = 0))

# Filter cell type pairs (Endo/DC)
df1 %>%
  filter(!cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC")) -> df1

#########################
# Network for each LR pair
#########################
## Define edges
dummy <- df1 %>%
  separate(cell_type_pair, sep = "->", into = c("from", "to")) %>%
  mutate(thickness = 0) %>%
  select(from, to, thickness) %>%
  distinct()


lrs <- unique(df1$LR)

df1_long <- df1


f <- function(Edges, pdfname, lr) {
  tryCatch(
    {
      pdf(pdfname)
      qgraph(Edges,
        esize = 5,
        theme = "colorblind", layout = "circle", # maximum=max_value,
        posCol = "blue",
        edge.labels = TRUE,
        cut = 7,
        title = sprintf("%s, number of patients", lr)
      )
      dev.off()
    },
    warning = function(w) {
      print(w)
    },
    error = function(e) {
      print(e)
    }
  )
}

for (lr in lrs) {
  # min_value <- min(c(df1sub$adjusted_mean_NormHRL, df1sub$adjusted_mean_NormHRH))
  # max_value <- max(c(df1sub$adjusted_mean_NormHRL, df1sub$adjusted_mean_NormHRH))

  df1_long %>%
    filter(LR == lr) %>%
    separate(cell_type_pair, sep = "->", into = c("from", "to")) %>%
    mutate(thickness = number_of_patient) %>%
    select(from, to, thickness) -> Edges

  bind_rows(Edges, dummy) %>%
    arrange(from) -> Edges

  pdfname <- file.path(
    path_outdir,
    sprintf("network_cell_type_%s_numberOfPatient.pdf", lr)
  )

  f(Edges, pdfname, lr)
}

# sessionInfo()
sessionInfo()
