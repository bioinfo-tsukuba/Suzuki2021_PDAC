################################################################################
# Initialization
################################################################################

options(repos = "https//cran.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, janitor)
system("mkdir -p results/Fig4/")

################################################################################
# Input and format
################################################################################

df_prognostic_lr <-
  read_csv("results/Fig1/LR_adjPval_meanHR_screened.csv", col_types = cols()) %>%
  mutate(prognosis = if_else(meanHR > 1, "poor", "good")) %>%
  select(LR, prognosis)

files <- system("find results/NATMI_each_patient/ExtractEdges/ | grep Edges_lrc2p.csv | sort", intern = TRUE)
df_natmi <- map2_dfr(files, seq_along(files), ~ read_csv(.x, col_types = cols()) %>% mutate(id = .y))

df_natmi_formatted <-
  df_natmi %>%
  clean_names() %>%
  select(id, 1:4, edge_average_expression_weight) %>%
  unite("CCI", c(sending_cluster, target_cluster), sep = "->") %>%
  unite("LR", c(ligand_symbol, receptor_symbol), sep = "->")

################################################################################
# Calculate threshold using Grade4 patients LR expression
################################################################################

df_natmi_grade4 <-
  df_natmi_formatted %>%
  # Grade4 patients
  filter(id == 1 | id == 3) %>%
  group_by(CCI, LR) %>%
  mutate(mean = mean(edge_average_expression_weight)) %>%
  select(!c(id, edge_average_expression_weight)) %>%
  distinct() %>%
  ungroup(LR) %>%
  mutate(median = median(mean)) %>%
  mutate(q1 = quantile(mean, 0.25)) %>%
  mutate(q3 = quantile(mean, 0.75))

################################################################################
# Join with df_prognostic_lr
################################################################################

df_result <-
  df_prognostic_lr %>%
  inner_join(df_natmi_grade4, by = "LR") %>%
  select(CCI, LR, prognosis, mean, median, q1, q3) %>%
  rename(mean_exp_weight_grade4 = mean,
    median_exp_weight_grade4 = median,
    q1_exp_weight_grade4 = q1,
    q3_exp_weight_grade4 = q3,
    )

################################################################################
# Output
################################################################################

write_csv(df_result, "results/Fig4/quantile_threshold_grade4.csv")
