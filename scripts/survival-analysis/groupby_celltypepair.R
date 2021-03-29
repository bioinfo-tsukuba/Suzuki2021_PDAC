# Report Peto & Peto modification of the Gehan-Wilcoxon test
system("mkdir -p results/MD3/")

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse, rvest, janitor)

df_survival <-
  read_csv("data/ICGC/survival_PAAD-US.csv.gz") %>%
  mutate(status = if_else(status == "alive", 0, 1)) # alive=0, dead=1

df_LR <-
  read_tsv("results/MD2/DiffEdges/tableForHeatmap.tsv") %>%
  # remove log2-transformed data because of Inf value
  filter(!str_detect(weight, "log2")) %>%
  filter(category != "All_edges_mean") %>%
  mutate(lr_pair = ligandreceptor_pair) %>%
  separate(ligandreceptor_pair, c("ligand", "receptor"), sep = "->") %>%
  pivot_longer(c(ligand, receptor), names_to = "lr", values_to = "gene")

report_pval <- function(data) {
  survdiff(Surv(time, status) ~ exp_bin, data = data, rho = 1) %>%
  glance() %>%
  pull(p.value)
}

#*+++++++++++++++++++++++++++++++++++++++++++++++++++
df_test <- df_LR %>% filter(category == "Appeared_mean", weight == "delta_edge_specificity_weight")
#*+++++++++++++++++++++++++++++++++++++++++++++++++++

df_pval <-
  df_survival %>%
  inner_join(df_test, key = "gene")

unique(df_pval$celltype_pair)
df_pval %>% select(celltype_pair, lr_pair) %>% filter(lr_pair == "A2M->LRP1") %>% count(celltype_pair)
tmp <- df_pval %>%
  group_by(celltype_pair, lr_pair) %>%
  # high and low by median value
  mutate(exp_sum = sum(exp)) %>%
  mutate(exp_bin = if_else(exp_sum > quantile(exp_sum, 0.5), "high", "low")) %>%
  # median survival time
  group_by(celltype_pair, exp_bin) %>%
  mutate(median_time = quantile(time, 0.5)) %>%
  # Report p-value
  nest(data = !c(celltype_pair))

tmp %>% filter(celltype_pair == "CAF->EMT") %>%
  unnest(data) %>% select(status, exp_bin) %>% distinct()
  mutate(pval = map_dbl(data, report_pval)) %>%
  filter(pval < 0.01) %>%
  unnest(data) %>%
  select(celltype_pair, pval, exp_bin, median_time) %>%
  distinct() %>%
  pivot_wider(c(celltype_pair, pval), names_from = "exp_bin", names_prefix = "median_day_", values_from = "median_time") %>%
  mutate(diff_day = median_day_high - median_day_low) %>%
  arrange(pval)

write_csv(df_pval, "results/MD3/survival_all_LR.csv")

