
library(tidyverse)
library(ggplot2)

#LR paris in each cell type
#All_edges_mean 

df_all_edge <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/MD2/DiffEdges/Delta_edges_lrc2p/All_edges_mean.csv")

df_all_edge %>% 
  group_by(`Sending cluster`, `Target cluster`) %>%
  summarise(n_lrpair = n()) -> df_all_edge_summary

write_csv(df_all_edge_summary, "~/Desktop/SSD/analysis/cell_type/df_all_edge_summary")

df_all_edge_summary %>%
  mutate(Cell_type_pair = paste0(`Sending cluster`, `Target cluster`)) %>%
  ggplot(aes(Cell_type_pair, n_lrpair))  +
  geom_bar(stat="identity") +
  theme(axis.text=element_text(size=3))-> g1

plot(g1)

ggsave("~/Desktop/SSD/analysis/cell_type/df_all_edge_summary.pdf", g1)

#Appeared
df_appeared <- read_csv('~/Desktop/SSD/results/MD2/DiffEdges/Delta_edges_lrc2p/Appeared_mean.csv')

df_appeared %>% 
  group_by(`Sending cluster`, `Target cluster`) %>%
  summarise(n_lrpair = n()) -> df_appeared_summary

write_csv(df_all_edge_summary, "~/Desktop/SSD/analysis/cell_type/df_appeared_summary")

df_appeared_summary %>%
  mutate(Cell_type_pair = paste0(`Sending cluster`, `Target cluster`)) %>%
  ggplot(aes(Cell_type_pair, n_lrpair))  +
  geom_bar(stat="identity") +
  theme(axis.text=element_text(size=3))-> g2

ggsave("~/Desktop/SSD/analysis/cell_type/df_appeared_summary.pdf", g2)

#Up-redurated
df_upregulated <- read_csv('~/Desktop/SSD/results/MD2/DiffEdges/Delta_edges_lrc2p/Up-regulated_mean.csv')

df_upregulated %>%
  group_by(`Sending cluster`, `Target cluster`) %>%
  summarise(n_lrpair = n()) -> df_upregulated_summary

write_csv(df_upregulated_summary, "~/Desktop/SSD/analysis/cell_type/df_upregulated_summary")

df_upregulated_summary %>%
  mutate(Cell_type_pair = paste0(`Sending cluster`, `Target cluster`)) %>%
  ggplot(aes(Cell_type_pair, n_lrpair))  +
  geom_bar(stat="identity") +
  theme(axis.text=element_text(size=3))-> g3

ggsave("~/Desktop/SSD/analysis/cell_type/df_upregulated_summary.pdf", g3)




