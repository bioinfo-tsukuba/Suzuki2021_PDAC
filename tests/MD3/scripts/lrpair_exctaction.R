if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, janitor, glue)

system("mkdir -p tests/MD3/data")

file_dir <- "results/MD2/DiffEdges/Delta_edges_lrc2p/Appeared_mean.csv"

df <-
  read_csv(file_dir, col_names = TRUE, col_types = cols()) %>%
  clean_names()

df_formatted <-
  df %>%
  mutate(celltype_pair = paste0(sending_cluster, "->", target_cluster)) %>%
  mutate(ligandreceptor_pair = paste0(ligand_symbol, "->", receptor_symbol)) %>%
  select(celltype_pair, ligandreceptor_pair,
    delta_edge_expression_weight, delta_edge_specificity_weight) %>%
  pivot_longer(-c(celltype_pair, ligandreceptor_pair),
    names_to = "weight", values_to = "value") %>%
    select(ligandreceptor_pair) %>%
    distinct()

write_csv(df_formatted, 'tests/MD3/data/lrpair.csv')