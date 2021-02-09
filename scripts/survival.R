# survminerパッケージを用いて生存曲線を描きます. P値も付けます. 

if (!require("BiocManager", quietly = TRUE)) install.packages("pacman")
pacman::p_load(survival, survminer, tidyverse)

df_raw <- read_csv("data/ICGC/survival.csv.gz")

# 適当にKRAS遺伝子をチョイスしました
df_gene <- df_raw %>% filter(gene == "KRAS")

# ggplot(df_gene, aes(x = gene, y = exp)) +
# geom_violin() +
# geom_boxplot() +
# geom_point()

exp_median <- quantile(df_gene$exp ,0.5)

df_plot <-
  df_gene %>%
  mutate(exp= if_else(exp > exp_median, "high", "low")) %>%
  mutate(status=if_else(status == "alive", 0, 1))

fit <- survfit(Surv(time, status) ~ exp, data = df_plot)

p <- ggsurvplot(fit, data = df_plot, pval = TRUE)
p