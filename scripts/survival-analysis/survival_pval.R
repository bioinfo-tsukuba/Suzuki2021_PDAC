# Report Peto & Peto modification of the Gehan-Wilcoxon test
system("mkdir -p results/MD3/")

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse)

df_survival <-
  read_csv("data/ICGC/survival.csv") %>%
  mutate(status = if_else(status == "alive", 0, 1)) # alive=0, dead=1

df_LR <-
  read_tsv("results/MD2/DiffEdges/tableForHeatmap.tsv")

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
    slice_max(delta_edge_specificity_weight, n = 20) %>%
    mutate(lr_pair = ligandreceptor_pair) %>%
    separate(ligandreceptor_pair, c("ligand","receptor")) %>%
    pivot_longer(-c(category, celltype_pair, lr_pair, delta_edge_specificity_weight),
      names_to = "ligandreceptor", values_to = "gene") %>%
    # Bind LR data and survival data
    inner_join(df_survival, by = "gene") %>%
    group_by(lr_pair) %>%
    # high and low by median value
    mutate(exp_bin = if_else(exp > quantile(exp, 0.5), "high", "low")) %>%
    nest(data = !c(category, lr_pair, delta_edge_specificity_weight)) %>%
    # Report p-value
    mutate(pval = map_dbl(data, report_pval)) %>%
    mutate(sig = if_else(pval < 0.05, TRUE, FALSE)) %>%
    select(category, lr_pair, delta_edge_specificity_weight, pval, sig)
  })

write_csv(df_pval, "results/MD3/survival_pval.csv")
