################################################################################
# Initialization
################################################################################

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse)
system("mkdir -p results/MD3/")

################################################################################
# Input and format
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  df_survival <- read_csv(args[1], col_types = cols())
  df_lr <- read_csv(args[2], col_names = "LR", col_types = cols())
} else {
  df_survival <- read_csv("data/ICGC/survival_PAAD-US.csv.gz", col_types = cols())
  df_lr <- read_csv("tests/MD3/data/LR.csv", col_names = "LR", col_types = cols())
}

df_survival <- df_survival %>%
  mutate(status = if_else(status == "alive", 0, 1)) # alive=0, dead=1

df_lr <- df_lr %>%
  mutate(lr_pair = LR) %>%
  separate(lr_pair, c("ligand", "receptor"), sep = "->") %>%
  pivot_longer(c(ligand, receptor), names_to = "lr", values_to = "gene")

################################################################################
# Function
################################################################################

report_pval <- function(data) {
  survdiff(Surv(time, status) ~ exp_bin, data = data, rho = 1) %>%
  glance() %>%
  pull(p.value)
}

report_hr <- function(data) {
  cox <- coxph(Surv(time, status) ~ exp_bin, data = data)
  exp(cox$coefficients)
}

################################################################################
# Main
################################################################################

df_pval <-
  inner_join(df_survival, df_lr, by = "gene") %>%
  group_by(LR, id) %>%
  # high and low by median value
  mutate(exp_sum = sum(exp)) %>%
  ungroup(id) %>%
  mutate(exp_bin = if_else(exp_sum > median(exp_sum), "high", "low")) %>%
  # median survival time
  group_by(LR, exp_bin) %>%
  mutate(median_time = median(time)) %>%
  select(LR, time, status, exp_bin, median_time) %>%
  distinct() %>%
  # Report P-value
  nest(data = !c(LR)) %>%
  mutate(pval = map_dbl(data, report_pval)) %>%
  # Report HR
  mutate(hr = map_dbl(data, report_hr)) %>%
  filter(pval < 0.05) %>%
  arrange(pval) %>%
  unnest(data) %>%
  select(LR, pval, hr, exp_bin, median_time) %>%
  distinct() %>%
  pivot_wider(c(LR, pval, hr), names_from = "exp_bin", names_prefix = "median_day_", values_from = "median_time") %>%
  mutate(diff_day = median_day_high - median_day_low) %>%
  arrange(pval)

################################################################################
# Output
################################################################################

cat(format_csv(df_pval))

