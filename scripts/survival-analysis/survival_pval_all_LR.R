# Report Peto & Peto modification of the Gehan-Wilcoxon test
system("mkdir -p results/MD3/")

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse, rvest, janitor)

df_survival <-
  read_csv("data/ICGC/survival.csv") %>%
  mutate(status = if_else(status == "alive", 0, 1)) # alive=0, dead=1

df_all_LR <-
  read_html("https://asrhou.github.io/NATMI/#2-connectomedb2020-") %>%
  html_table() %>%
  flatten_dfr() %>%
  clean_names() %>%
  mutate(pair = str_c(ligand_gene_symbol, receptor_gene_symbol, sep = "->")) %>%
  select(pair, ligand_gene_symbol, receptor_gene_symbol) %>%
  pivot_longer(-pair, names_to="LR", values_to="gene")

report_pval <- function(data) {
  survdiff(Surv(time, status) ~ exp_bin, data = data, rho = 1) %>%
  glance() %>%
  pull(p.value)
}

df_pval <-
  df_survival %>%
  inner_join(df_all_LR, key = "gene") %>%
  group_by(pair, id, status, time) %>%
  # high and low by median value
  summarize(exp_sum = sum(exp)) %>%
  ungroup(id, status, time) %>%
  mutate(exp_bin = if_else(exp_sum > quantile(exp_sum, 0.5), "high", "low")) %>%
  # median survival time
  group_by(pair, exp_bin) %>%
  mutate(median_time = quantile(time, 0.5)) %>%
  # Report p-value
  nest(data = !c(pair)) %>%
  mutate(pval = map_dbl(data, report_pval)) %>%
  filter(pval < 0.01) %>%
  unnest(data) %>%
  select(pair, pval, exp_bin, median_time) %>%
  distinct() %>%
  pivot_wider(c(pair, pval), names_from = "exp_bin", names_prefix = "median_day_", values_from = "median_time") %>%
  mutate(diff_day = median_day_high - median_day_low) %>%
  arrange(pval)

write_csv(df_pval, "results/MD3/survival_all_LR.csv")

