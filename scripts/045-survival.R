# survival.shで整形したデータを用いて生存曲線を描きます. P値も付けます.

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, tidyverse)

df_raw <- read_csv("data/ICGC/survival.csv")

# 遺伝子をチョイスします（任意）
df_gene <- df_raw %>% filter(gene == "WNT4")
df_gene <- df_raw %>% filter(gene == "CD274") # PDL1
df_gene <- df_raw %>% filter(gene == "PDCD1") # PD1

# ggplot(df_gene, aes(x = gene, y = exp)) +
# geom_violin() +
# geom_boxplot() +
# geom_point()

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
