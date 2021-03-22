if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, janitor, glue)

file_dir <- "results/MD2/DiffEdges/Delta_edges_lrc2p/"
files <- dir(file_dir, pattern = "*.csv")
files <- files[!str_detect(files, "Stable")]

df <- map_dfr(files, function(.x) {
  read_csv(paste0(file_dir, .x), col_names = TRUE, col_types = cols()) %>%
  mutate(category = str_remove(.x, ".csv")) %>%
  clean_names()
  })

df_formatted <- df %>%
  mutate(celltype_pair = paste0(sending_cluster, "->", target_cluster)) %>%
  mutate(ligandreceptor_pair = paste0(ligand_symbol, "->", receptor_symbol)) %>%
  select(category, celltype_pair, ligandreceptor_pair,
    delta_edge_expression_weight, delta_edge_specificity_weight,
    log2_transformed_fold_change_of_edge_expression_weight,
    log2_transformed_fold_change_of_edge_specificity_weight) %>%
  pivot_longer(-c(category, celltype_pair, ligandreceptor_pair),
    names_to = "weight", values_to = "value")

write_tsv(df_formatted, 'results/MD2/DiffEdges/tableForHeatmap.tsv')

walk(df_formatted$weight %>% unique, ~ {
  df_formatted %>%
    filter(weight == .x) %>%
    group_by(category) %>%
    slice_max(value, n = 20) %>%
    ggplot(aes(x = celltype_pair, y = ligandreceptor_pair)) +
    geom_raster(aes(fill = value)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_classic() +
    theme(text = element_text(size = 12)) +
    labs(title = "Top 20 pairs",
      fill = glue({.x}),
      x = NULL,
      y = NULL) +
    facet_wrap(~ category, scales = "free")
    ggsave(glue("results/MD2/DiffEdges/heatmap_top20_{.x}.pdf"), width = 30)
})
