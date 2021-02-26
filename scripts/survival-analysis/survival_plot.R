# Survival plot
system("mkdir -p results/MD3/")

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, patchwork, survival, survminer)

df_survival <-
  read_csv("data/ICGC/survival.csv") %>%
  mutate(status = if_else(status == "alive", 0, 1)) # alive=0, dead=1

df_LR <-
  read_tsv("results/MD2/DiffEdges/tableForHeatmap.tsv")

df_plot <-
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
    mutate(exp_bin = if_else(exp > quantile(exp, 0.5), "high", "low"))
  }) %>%
  as.data.frame()

walk(df_LR$category %>% unique, function(.x) {
  df_tmp <- df_plot %>% filter(category == .x)
  fit <- survfit(Surv(time, status) ~ exp_bin, data = df_tmp)
  g <- ggsurvplot_facet(fit, df_tmp, facet.by = "lr_pair", pval = TRUE)
  ggsave(sprintf("results/MD3/survival_plot_%s.pdf", .x), g, dpi = 300, width = 30, height = 20, limitsize = FALSE)
})