# survival.shで整形したデータを用いて生存曲線を描きます. P値も付けます.

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, tidyverse)

df_raw <- read_csv("data/ICGC/survival.csv")

# 任意の遺伝子をチョイスします
df_gene <- df_raw %>% filter(gene == "CD274") # PDL1
df_gene <- df_raw %>% filter(gene == "PDCD1") # PD1
df_gene <- df_raw %>% filter(gene == "CXCL12")
df_gene <- df_raw %>% filter(gene == "CXCR4")


# ggplot(df_gene, aes(x = gene, y = exp)) +
#   geom_violin() +
#   geom_boxplot() +
#   geom_point()


# Genome medicine[https://doi.org/10.1186/s13073-020-00776-9]の論文に倣って,
# 中央値をhighとlowの閾値とします

exp_median <- quantile(df_gene$exp, 0.5)

df_plot <-
  df_gene %>%
  mutate(exp= if_else(exp > exp_median, "high", "low")) %>%
  mutate(status=if_else(status == "alive", 0, 1))

fit <- survfit(Surv(time, status) ~ exp, data = df_plot)

p <- ggsurvplot(fit, data = df_plot, pval = TRUE)
p


### facetを試そうとしましたが, エラーです…

# genes <- c("CD274","PDCD1","CXCL12","CXCR4")
# df_gene <- df_raw %>% filter(gene %in% genes)

# Genome medicine[https://doi.org/10.1186/s13073-020-00776-9]の論文に倣って,
# 中央値をhighとlowの閾値とします

# exp_median <- quantile(df_gene$exp, 0.5)

# df_plot <-
#   df_gene %>%
#   mutate(status = if_else(status == "alive", 0, 1)) %>% # alive=0, dead=1
#   group_by(gene) %>%
#   mutate(exp = if_else(exp < quantile(exp, 0.5), 0, 1)) %>% # low=0, high=1
#   mutate(gene = as.factor(gene)) %>%
#   ungroup(gene)

# fit <- survfit(Surv(time, status) ~ sex, data = as_tibble(colon))
# as_tibble(colon) %>% select(time, status, sex, rx)
# ggsurvplot_facet(fit, colon, facet.by = "rx", pval = TRUE)

# fit1 <- survfit(Surv(time, status) ~ exp, data = df_plot)
# df_plot %>% select(time, status, exp, gene)
# ggsurvplot_facet(fit1, df_plot, facet.by = "gene", pval = TRUE)

# tmp <- head(df_plot)

# fit1 <- survfit(Surv(time, status) ~ exp, data = tmp)
# tmp %>% select(time, status, exp, gene)
# ggsurvplot_facet(fit1, tmp, facet.by = "gene", pval = TRUE)

#====================================
# > ggsurvplot_facet(fit1, tmp, facet.by = "gene", pval = TRUE)
#  if (xi > xj) 1L else -1L でエラー: 
#    TRUE/FALSE が必要なところが欠損値です 
#  追加情報:  警告メッセージ: 
#  Ops.factor(xi, xj) で:  ‘>’ not meaningful for factors
# > 
# > tmp$gene
# [1] CD274  CXCL12 CXCR4  PDCD1  CD274  CXCL12
# Levels: CD274 CXCL12 CXCR4 PDCD1