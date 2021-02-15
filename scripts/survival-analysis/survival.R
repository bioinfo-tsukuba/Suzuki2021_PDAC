# Survival plot

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, tidyverse)

df_raw <- read_csv("data/ICGC/survival.csv")

# Target gene lists

genes <- c("CD274","PDCD1","CXCL12","CXCR4")
df_gene <- df_raw %>% filter(gene %in% genes)

# high and low by median value

df_plot <-
  df_gene %>%
  mutate(status = if_else(status == "alive", 0, 1)) %>% # alive=0, dead=1
  mutate(gene = as.factor(gene)) %>%
  group_by(gene) %>%
  mutate(exp_bin = if_else(exp > quantile(exp, 0.5), "high", "low")) %>%
  as.data.frame()

fit <- survfit(Surv(time, status) ~ exp_bin, data = df_plot)
ggsurvplot_facet(fit, df_plot, facet.by = "gene", pval = TRUE)

ggsave("analysis/survival/testplot.png", dpi = 600)

## FYI: Gene expression
# ggplot(df_plot, aes(x = "", y = exp)) +
#   geom_violin() +
#   geom_point() +
#   facet_wrap("gene", scales = "free")
