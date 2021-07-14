################################################################################
# Initialization
################################################################################

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, svglite)

################################################################################
# Input and format
################################################################################

df <- read_csv("results/Fig1/LR_HR_adjPval.csv")

df_plot <-
  df %>%
  select(!starts_with(c("Pval", "meta_Pval"))) %>%
  pivot_longer(!c(LR, adjPval), names_to = "group", values_to = "HR") %>%
  mutate(mlogPval = -log10(adjPval)) %>%
  group_by(LR) %>%
  mutate(mean_HR = mean(HR)) %>%
  select(LR, mlogPval, mean_HR) %>%
  distinct() %>%
  mutate(label = case_when(
    mlogPval > 1 && mean_HR > 1 ~ "poor",
    mlogPval > 1 && mean_HR < 1 ~ "good",
    TRUE ~ "non"
  )) %>%
  ungroup()

df_plot %>%
  group_by(label) %>%
  count(label)

p <- ggplot(df_plot, aes(x = mean_HR, y = mlogPval, color = label)) +
  geom_point() +
  scale_color_manual(values = c("#ff5050", "gray", "#0066ff")) +
  labs(x = "mean HR", y = "-log10 P-value") +
  theme(text = element_text(size = 16, family = "Arial")) +
  theme_bw()

ggsave("results/Fig1/dotplot.svg", p, width = 5, height = 3)
