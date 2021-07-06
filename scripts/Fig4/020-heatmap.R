################################################################################
# Initialization
################################################################################

options(repos = "https//cran.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, janitor, patchwork)
system("mkdir -p results/Fig4/")

################################################################################
# Input and format
################################################################################

df_lr <- read_csv("results/Fig4/quantile_threshold_grade4.csv", col_types = cols())

files <- system("find results/NATMI_each_patient/ExtractEdges/ | grep Edges_lrc2p.csv | sort", intern = TRUE)
patients_grade <-
  read_csv("data/patients_info.csv", col_types = cols()) %>%
  filter(Primary_Metastasis == "Primary") %>%
  pull(Grade)
df_natmi <- map2_dfr(files, patients_grade, ~ read_csv(.x, col_types = cols()) %>% mutate(grade = .y))

df_natmi_formatted <-
  df_natmi %>%
  clean_names() %>%
  select(grade, 1:4, edge_average_expression_weight) %>%
  unite("CCI", c(sending_cluster, target_cluster), sep = "->") %>%
  unite("LR", c(ligand_symbol, receptor_symbol), sep = "->") %>%
  group_by(grade, CCI, LR) %>%
  summarize(mean_weight = mean(edge_average_expression_weight)) %>%
  ungroup() %>%
  # ENDOとDCの組み合わせを除きます.
  filter(!CCI %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC"))


################################################################################
# Filterしない場合のヒートマップを描出します
################################################################################

df_plot <-
  df_natmi_formatted %>%
  inner_join(df_lr, by = c("CCI", "LR")) %>%
  filter(prognosis == "poor")

df_plot %>%
  select(LR) %>%
  n_distinct() # 55

plot_tile <- function(data, title) {
  ggplot(data, aes(x = grade, y = fct_rev(LR), fill = mean_weight)) +
    geom_tile() +
    scale_fill_gradient(limits = c(0, max(data$mean_weight)), low = "#EEEEEE", high = "blue") +
    labs(title = title, x = "Grade", y = "LR pair", fill = "Score")
}

g_no_filter <-
  df_plot %>%
  group_nest(CCI) %>%
  mutate(g = map2(data, CCI, plot_tile)) %>%
  pull(g) %>%
  wrap_plots()

ggsave("results/Fig4/heatmap_no_filter.pdf", g_no_filter, width = 30, height = 30)

################################################################################
# Grade4での平均の発現が中央値よりも高いペアを抽出します
################################################################################

df_plot <-
  df_natmi_formatted %>%
  inner_join(df_lr, by = c("CCI", "LR")) %>%
  filter(mean_exp_weight_grade4 > median_exp_weight_grade4) %>%
  filter(prognosis == "poor")

df_plot %>%
  select(LR) %>%
  n_distinct() # 47

g_grade4 <-
  df_plot %>%
  group_nest(CCI) %>%
  mutate(g = map2(data, CCI, plot_tile)) %>%
  pull(g) %>%
  wrap_plots()

ggsave("results/Fig4/heatmap_median_grade4.pdf", g_grade4, width = 30, height = 30)


################################################################################
# Grade4での平均の発現が中央値よりも高いペアを抽出します
# 2<3<4と発現が上昇しているLRペアを抽出します.
################################################################################

df_increase <-
  df_natmi_formatted %>%
  group_by(CCI, LR) %>%
  mutate(rank = rank(mean_weight) + 1) %>%
  count(trend = (grade == rank)) %>%
  filter(trend == TRUE, n == 3) %>%
  select(CCI, LR)

df_plot <-
  df_natmi_formatted %>%
  inner_join(df_lr, by = c("CCI", "LR")) %>%
  inner_join(df_increase, by = c("CCI", "LR")) %>%
  filter(mean_exp_weight_grade4 > median_exp_weight_grade4) %>%
  filter(prognosis == "poor")

df_plot %>%
  select(LR) %>%
  n_distinct() # 23

g <-
  df_plot %>%
  group_nest(CCI) %>%
  mutate(g = map2(data, CCI, plot_tile)) %>%
  pull(g) %>%
  wrap_plots()

ggsave("results/Fig4/heatmap_median_grade4_trend.pdf", g, width = 30, height = 30)


################################################################################
# Grade4での平均の発現が75%値よりも高いペアを抽出します
################################################################################

df_plot <-
  df_natmi_formatted %>%
  inner_join(df_lr, by = c("CCI", "LR")) %>%
  filter(mean_exp_weight_grade4 > q3_exp_weight_grade4) %>%
  filter(prognosis == "poor")

df_plot %>%
  select(LR) %>%
  n_distinct() # 33

g_grade4 <-
  df_plot %>%
  group_nest(CCI) %>%
  mutate(g = map2(data, CCI, plot_tile)) %>%
  pull(g) %>%
  wrap_plots()

ggsave("results/Fig4/heatmap_q3_grade4.pdf", g_grade4, width = 30, height = 30)
