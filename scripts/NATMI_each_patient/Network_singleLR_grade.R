
path_result <- "results/NATMI_each_patient/Dataframe_AtoE/Datafram_filtered_grade.csv"
path_metadata <- "data/WeiLin_pdac10/Patient_data.csv"
path_outdir <- "results/NATMI_each_patient/Network_singleLR_grade"

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

# Read metadata
dfmeta <- read_csv(path_metadata) %>%
  filter(Primary_or_Metasitasis == "Primary")


#########################
# Network for each grade
#########################
## Define edges
grades <- sort(unique(dfmeta$Grade))

dummy <- df1 %>%
  separate(cell_type_pair, sep = "->", into = c("from", "to")) %>%
  mutate(thickness = 0) %>%
  select(from, to, thickness) %>%
  distinct()


lrs <- unique(df1$LR)

df1 %>%
  pivot_longer(
    c(`2`, `3`, `4`),
    names_to = "Grade",
    values_to = "Averaged_mean_expression_weight"
  ) -> df1_long


f <- function(Edges, pdfname, lr, grade) {
  tryCatch(
    {
      pdf(pdfname)
      qgraph(Edges,
        esize = 5,
        theme = "colorblind", layout = "circle", # maximum=max_value,
        title = sprintf("%s, Grade %s", lr, grade)
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
  for (grade in grades) {

    # min_value <- min(c(df1sub$adjusted_mean_NormHRL, df1sub$adjusted_mean_NormHRH))
    # max_value <- max(c(df1sub$adjusted_mean_NormHRL, df1sub$adjusted_mean_NormHRH))

    # Define edges, HR<1
    df1_long %>%
      filter(LR == lr & Grade == grade) %>%
      separate(cell_type_pair, sep = "->", into = c("from", "to")) %>%
      mutate(thickness = Averaged_mean_expression_weight) %>%
      select(from, to, thickness) -> Edges

    bind_rows(Edges, dummy) %>%
      arrange(from) -> Edges

    pdfname <- file.path(
      path_outdir,
      sprintf("network_cell_type_%s_Grade%s.pdf", lr, grade)
    )

    f(Edges, pdfname, lr, grade)
  }
}

# sessionInfo()
sessionInfo()
