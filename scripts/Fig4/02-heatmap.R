################################################################################
# Initialization
################################################################################
options(repos = "https//cran.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, janitor, gridExtra)

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
  filter(CCI %in% c("ETC->TAM", "CAF->ETC", "ETC->CAF"))

plot_tile <- function(data, title) {
  ggplot(data, aes(x = grade, y = fct_rev(LR), fill = mean_weight)) +
    geom_tile() +
    scale_fill_gradient(limits = c(0, max(data$mean_weight)), low = "#EEEEEE", high = "blue") +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_text(size = 8, family = "Helvetica", color = "black"),
      axis.text = element_text(size = 6, family = "Helvetica", color = "black"),
      legend.text = element_text(size = 6, family = "Helvetica", color = "black"),
      legend.title = element_text(size = 6, family = "Helvetica", color = "black"),
      strip.text = element_text(size = 6, family = "Helvetica", color = "black"),
      legend.key.size = unit(0.05, "inch"),
      panel.spacing.x = unit(0, "lines"),
      panel.grid.minor = element_blank()
    ) +
    labs(title = title, y = "", x = "", fill = "Score")
}

################################################################################
# Main
################################################################################

df_plot <-
  df_natmi_formatted %>%
  inner_join(df_lr, by = c("CCI", "LR")) %>%
  filter(mean_exp_weight_grade4 > q3_exp_weight_grade4) %>%
  filter(prognosis == "poor")

g_grade4 <-
  df_plot %>%
  select(CCI, LR, grade, mean_weight) %>%
  complete(grade, nesting(CCI, LR), fill = list(mean_weight = 0)) %>%
  group_nest(CCI) %>%
  mutate(g = map2(data, CCI, plot_tile)) %>%
  pull(g)

g <- marrangeGrob(g_grade4, nrow = 1, ncol = 3)

ggsave("results/Fig4/heatmap_q3_grade4.pdf", g, width = 174, height = 50, units = "mm")