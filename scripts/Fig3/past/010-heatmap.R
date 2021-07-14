################################################################################
# Initialization
################################################################################

options(repos = "https//cran.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, janitor, ggh4x, ggdendro)
system("mkdir -p results/Fig3/")

################################################################################
# Input and format
################################################################################

df_prognostic_lr <-
  read_csv("results/Fig1/LR_adjPval_meanHR_screened.csv", col_types = cols()) %>%
  mutate(prognosis = if_else(meanHR > 1, "poor", "good")) %>%
  select(LR, prognosis)

files <- system("find results/NATMI_each_patient/ExtractEdges/ | grep Edges_lrc2p.csv | sort", intern = TRUE)
df_natmi <- map2_dfr(files, seq_along(files), ~ read_csv(.x, col_types = cols()) %>% mutate(id = .y))

df_natmi_formatted <-
  df_natmi %>%
  clean_names() %>%
  select(id, 1:4, edge_average_expression_weight) %>%
  unite("CCI", c(sending_cluster, target_cluster), sep = "->") %>%
  unite("LR", c(ligand_symbol, receptor_symbol), sep = "->")

################################################################################
# Count patients number
################################################################################

df_patients_num <-
  inner_join(df_prognostic_lr, df_natmi_formatted, key = "LR") %>%
  filter(prognosis == "poor") %>%
  group_by(LR, CCI) %>%
  mutate(n = length(id)) %>%
  select(LR, CCI, n) %>%
  ungroup() %>%
  distinct()

plot_heatmap <- function(data) {
  tmp_wider <-
    data %>%
    pivot_wider(names_from = c(CCI), values_from = n, values_fill = 0L)

  mat_n <- as.matrix(tmp_wider[, -1])
  colnames(mat_n) <- colnames(tmp_wider[, -1])
  rownames(mat_n) <- pull(tmp_wider, LR)

  clust_LR <- hclust(dist(mat_n), "ward.D2")
  clust_CCI <- hclust(dist(t(mat_n)), "ward.D2")

  ggplot(data, aes(x = CCI, y = LR, fill = n)) +
    geom_tile() +
    scale_x_dendrogram(hclust = clust_CCI) +
    scale_y_dendrogram(hclust = clust_LR) +
    scale_fill_gradient(limits = c(0, 10)) +
    theme_bw(base_size = 8)
}

################################################################################
# Plot heatmap
################################################################################

g <- plot_heatmap(df_patients_num)
ggsave("results/Fig3/heatmap_patients_num.pdf", g, width = 30, height = 15)



df_patients_num_filter <-
  df_patients_num %>%
  filter(n > 7)

g_fiter <- plot_heatmap(df_patients_num_filter)
ggsave("results/Fig3/heatmap_patients_8.pdf", g_fiter, width = 30, height = 15)
