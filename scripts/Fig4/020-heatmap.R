################################################################################
# Initialization
################################################################################

options(repos= "https//cran.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, janitor, patchwork)
system("mkdir -p results/Fig4/")

################################################################################
# Input and format
################################################################################

df_lr <- read_csv("results/Fig4/median_threshold_grade4.csv", col_types = cols())

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
  ungroup()

################################################################################
# Filter LR:
# (1) ENDOとDCの組み合わせを除きます.
# (2) 2<3<4と発現が上昇しているLRペアを抽出します.
################################################################################

df_filtered <-
  df_natmi_formatted %>%
  filter(!CCI %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC")) %>%
  group_by(CCI, LR) %>%
  mutate(rank = rank(mean_weight) + 1) %>%
  count(trend = (grade == rank)) %>%
  filter(trend == TRUE, n == 3) %>%
  select(CCI, LR)

################################################################################
# Join dataframes
# (1) 予後に関わるLRペア群df_lrと結合します
# (2) 前項のdf_filteredと結合します
# (3) Grade4での平均の発現が中央値よりも高いペアを抽出します
# (4) 予後不良因子のみを抽出します
################################################################################

df_plot <-
  df_natmi_formatted %>%
  inner_join(df_lr, by = c("CCI", "LR")) %>%
  inner_join(df_filtered, by = c("CCI", "LR")) %>%
  filter(candidate == TRUE) %>%
  filter(prognosis == "poor")

################################################################################
# Plot
################################################################################

plot_tile <- function(data) {
  ggplot(data, aes(x = grade, y = fct_rev(LR), fill = mean_weight)) +
    geom_tile() +
    scale_fill_gradient(limits=c(0, max(data$mean_weight))) +
    labs(title = .x, x = "Grade", y = "LR pair", fill = "Score")
}

g <-
  df_plot %>%
  group_nest(CCI) %>%
  mutate(g = map(data, plot_tile)) %>%
  pull(g) %>%
  wrap_plots()

ggsave("results/Fig4/heatmap.pdf", g, width = 30, height = 30)
