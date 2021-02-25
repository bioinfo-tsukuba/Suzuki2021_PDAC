# Survival plot
system("mkdir -p results/MD3/")

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, broom, tidyverse)

df_survival <- read_csv("data/ICGC/survival.csv") %>%
   mutate(status = if_else(status == "alive", 0, 1))  # alive=0, dead=1
df_LR <-
  read_tsv("results/MD2/DiffEdges/tableForHeatmap.tsv")

df_LR_top <-
  df_LR %>%
  group_by(category) %>%
  slice_max(delta_edge_specificity_weight, n = 20)

df_nested <-
  df_LR_top %>%
  mutate(lr_pair = ligandreceptor_pair) %>%
  separate(ligandreceptor_pair, c("ligand","receptor")) %>%
  pivot_longer(-c(category, celltype_pair, lr_pair, delta_edge_specificity_weight),
    names_to = "ligandreceptor", values_to = "gene") %>%
  nest_join(df_survival)

# high and low by median value

df_plot <-
  df_nested %>%
  rowwise() %>%

  group_by(LRpair) %>%
  mutate(exp_bin = if_else(exp > quantile(exp, 0.5), "high", "low")) %>%
  as.data.frame()

p_val <- survdiff(Surv(time, status) ~ exp_bin, data = df_plot, rho = 1) %>%
  glance() %>%
  pull(p.value)


fit <- survfit(Surv(time, status) ~ exp_bin, data = df_plot)
g <- ggsurvplot_facet(fit, df_plot, facet.by = "LRpair", pval = TRUE)

ggsave("results/MD3/survival_curve.pdf", g, dpi = 300, width = 20, height = 20)

## FYI: Gene expression
# ggplot(df_plot, aes(x = "", y = exp)) +
#   geom_violin() +
#   geom_point() +
#   facet_wrap("gene", scales = "free")
