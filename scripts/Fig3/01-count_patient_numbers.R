################################################################################
# Initialization
################################################################################
options(repos = "https//cran.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse)
system("mkdir -p results/Fig3/")

################################################################################
# Import and format
################################################################################

df_survival <- read_csv("results/survival_meta_analysis/NATMI_meta_filtered_Storey.csv")
df_survival_LR <- select(.data = df_survival, LR)

files <- system("find results/NATMI_each_patient/ExtractEdges/ | grep Edges_lrc2p.csv | sort", intern = TRUE)
patients <- 1:10 %>%
  str_pad(2, pad = c("0")) %>%
  str_c("P", .)

df_result <- map2_dfr(files, patients, ~
read_csv(.x, col_types = cols()) %>%
  mutate(Patient = .y) %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") %>%
  left_join(df_survival_LR, ., by = "LR"))

################################################################################
# Main
################################################################################

df_summary <-
  df_result %>%
  filter(!is.na(Patient)) %>%
  group_by(LR, cell_type_pair) %>%
  summarize(number_of_patient = n())

write_csv(df_summary, "results/Fig3/summary_number_of_patient.csv")
