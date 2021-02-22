# Survival plot

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, tidyverse)

df_raw <- read_csv("data/ICGC/survival.csv")
df_LR <-
  read_tsv("results/MD2/DiffEdges/upregulatedLR.tsv") %>%
  select(`Cell-type pair`, `Ligand-receptor pair`)

# Target gene lists
genes <- c("SST->SSTR2","C1QTNF1->AVPR2")

df_genes <-
  map_dfr(genes, function(x) {
    df_raw %>%
    filter(gene %in% str_split(x, "->", simplify=TRUE)) %>%
    mutate(LRpair = as.factor(x))}
    )

# high and low by median value

df_plot <-
  df_genes %>%
  mutate(status = if_else(status == "alive", 0, 1)) %>% # alive=0, dead=1
  group_by(LRpair) %>%
  mutate(exp_bin = if_else(exp > quantile(exp, 0.5), "high", "low")) %>%
  as.data.frame()

fit <- survfit(Surv(time, status) ~ exp_bin, data = df_plot)
ggsurvplot_facet(fit, df_plot, facet.by = "LRpair", pval = TRUE)

ggsave("analysis/survival/testplot.png", dpi = 600)

## FYI: Gene expression
# ggplot(df_plot, aes(x = "", y = exp)) +
#   geom_violin() +
#   geom_point() +
#   facet_wrap("gene", scales = "free")
