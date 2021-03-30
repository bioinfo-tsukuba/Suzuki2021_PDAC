# Report Peto & Peto modification of the Gehan-Wilcoxon test
system("mkdir -p results/MD3/")

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse)

df_survival <-
  read_csv("data/ICGC/survival.csv") %>%
  mutate(status = if_else(status == "alive", 0, 1)) # alive=0, dead=1

df_LR <-
  read_tsv("results/MD2/DiffEdges/tableForHeatmap.tsv") %>%
  filter(!str_detect(weight, "log2")) # remove log2-transformed data because of Inf value

report_pval <- function(data) {
  survdiff(Surv(time, status) ~ exp_bin, data = data, rho = 1) %>%
  glance() %>%
  pull(p.value)
}

df_pval <-
  map_dfr(df_LR$category %>% unique, function(.x) {
    df_LR %>%
    filter(category == .x) %>%
    # Extract top20
    group_by(weight) %>%
    slice_max(value, n = 20) %>%
    mutate(lr_pair = ligandreceptor_pair) %>%
    separate(ligandreceptor_pair, c("ligand","receptor")) %>%
    pivot_longer(-c(category, celltype_pair, lr_pair, weight, value),
      names_to = "ligandreceptor", values_to = "gene") %>%
    # Bind LR data and survival data
    inner_join(df_survival, by = "gene") %>%
    group_by(lr_pair, id) %>%
    # high and low by median value
    mutate(exp_sum = sum(exp)) %>%
    ungroup(id) %>%
    mutate(exp_bin = if_else(exp_sum > quantile(exp_sum, 0.5), "high", "low")) %>%
    # Report p-value
    nest(data = !c(category, celltype_pair, lr_pair, weight, value)) %>%
    mutate(pval = map_dbl(data, report_pval)) %>%
    filter(pval < 0.05) %>%
    select(category, celltype_pair, weight, lr_pair, pval)
  }) %>%
  arrange(pval)

write_csv(df_pval, "results/MD3/survival_pval_top20.csv")

