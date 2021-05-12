# Report Peto & Peto modification of the Gehan-Wilcoxon test
system("mkdir -p results/MD3/")

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse)

################################################################################
# Input
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  df_survival <- read_csv(args[1])
  df_lr <- read_csv(args[2])
} else {
  df_survival <- read_csv("data/ICGC/survival_PAAD-US.csv.gz", col_types = cols())
  df_lr <- read_csv("tests/MD3/data/lrpair.csv", col_types = cols())
}

df_survival <- df_survival %>%
  mutate(status = if_else(status == "alive", 0, 1)) # alive=0, dead=1

df_lr <- df_lr %>%
  mutate(lr_pair = ligandreceptor_pair) %>%
  separate(ligandreceptor_pair, c("ligand", "receptor"), sep = "->") %>%
  pivot_longer(c(ligand, receptor), names_to = "lr", values_to = "gene")

report_pval <- function(data) {
  survdiff(Surv(time, status) ~ exp_bin, data = data, rho = 1) %>%
    glance() %>%
    pull(p.value)
}

df_pval <-
  inner_join(df_survival, df_lr, by = "gene") %>%
  group_by(lr_pair, id) %>%
  # high and low by median value
  mutate(exp_sum = sum(exp)) %>%
  ungroup(id) %>%
  mutate(exp_bin = if_else(exp_sum > quantile(exp_sum, 0.5), "high", "low")) %>%
  # median survival time
  group_by(lr_pair, exp_bin) %>%
  mutate(median_time = quantile(time, 0.5)) %>%
  select(lr_pair, time, status, exp_bin, median_time) %>%
  distinct() %>%
  # Report p-value
  nest(data = !c(lr_pair)) %>%
  mutate(pval = map_dbl(data, report_pval)) %>%
  filter(pval < 0.05) %>%
  arrange(pval) %>%
  unnest(data) %>%
  select(lr_pair, pval, exp_bin, median_time) %>%
  distinct() %>%
  pivot_wider(c(lr_pair, pval), names_from = "exp_bin", names_prefix = "median_day_", values_from = "median_time") %>%
  mutate(diff_day = median_day_high - median_day_low) %>%
  arrange(pval)

cat(format_csv(df_pval))
# write_csv(df_pval, "results/MD3/survival_pval_top20.csv")
