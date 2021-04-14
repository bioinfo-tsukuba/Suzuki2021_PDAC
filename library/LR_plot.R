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
  df_lr <- read_csv(args[2], col_types = cols())
} else {
  df_survival <- read_csv("data/ICGC/survival_PAAD-US.csv.gz", col_types = cols())
  df_lr <- read_csv("tests/MD3/data/LR_top10.csv", col_types = cols())
}

df_survival <- df_survival %>%
  mutate(status = if_else(status == "alive", 0, 1))# alive=0, dead=1

df_lr <- df_lr %>%
  mutate(lr_pair = LR) %>%
  separate(lr_pair, c("ligand", "receptor"), sep = "->") %>%
  pivot_longer(c(ligand, receptor), names_to = "lr", values_to = "gene")

################################################################################
# Main
################################################################################

df_plot <-
  inner_join(df_survival, df_lr, by = "gene") %>%
  group_by(LR, id) %>%
  # high and low by median value
  mutate(exp_sum = sum(exp)) %>%
  ungroup(id) %>%
  mutate(exp_bin = if_else(exp_sum > median(exp_sum), "high", "low")) %>%
  select(LR, time, status, exp_bin) %>%
  distinct() %>%
  as.data.frame()

fit <- survfit(Surv(time, status) ~ exp_bin, data = df_plot)
g <- ggsurvplot_facet(fit, df_plot, facet.by = "LR", pval = TRUE)

ggsave("results/MD3/LR.pdf", g, dpi = 300, width = 20, height = 20)
