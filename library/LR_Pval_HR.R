# Calculate P-val and Hazard ratio

################################################################################
# Initialization
################################################################################

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse)

################################################################################
# Input and format
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  df_survival <- read_csv(args[1], col_types = cols())
  df_lr <- read_csv(args[2], col_names = "LR", col_types = cols())
} else {
  df_survival <- read_csv("data/ICGC/survival_PAAD-US.csv.gz", col_types = cols())
  df_lr <- read_csv("data/NATMI_LR.csv", col_names = "LR", col_types = cols())
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
  1 / exp(cox$coefficients)
}

################################################################################
# Main
################################################################################

df_pval <-
  inner_join(df_survival, df_lr, by = "gene") %>%
  select(id, status, time, gene, exp, LR) %>%
  group_by(id, LR) %>%
  # high and low by median value
  mutate(exp_sum = sum(exp)) %>%
  ungroup(id) %>%
  mutate(exp_bin = if_else(exp_sum > median(exp_sum), "high", "low")) %>%
  select(LR, time, status, exp_bin) %>%
  distinct() %>%
  # Report P-value
  nest(data = !c(LR)) %>%
  mutate(pval = map_dbl(data, report_pval)) %>%
  # Report HR
  mutate(hr = map_dbl(data, report_hr)) %>%
  arrange(pval) %>%
  select(LR, pval, hr) %>%
  distinct() %>%
  arrange(pval)

################################################################################
# Output
################################################################################

cat(format_csv(df_pval))
