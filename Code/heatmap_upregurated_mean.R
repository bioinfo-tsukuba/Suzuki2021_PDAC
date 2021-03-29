#heatmap for UP-regulated_mean.csv
library(tidyverse)
path_upregulated <- "results/MD2/DiffEdges/Delta_edges_lrc2p/UP-regulated_mean.csv"
df1 <- read_csv(path_upregulated, col_names=TRUE)

df1 %>% 
  mutate(`Cell-type pair`=paste0(`Sending cluster`, "->", `Target cluster`)) %>%
  mutate(`Ligand-receptor pair`=paste0(`Ligand symbol`, "->", `Receptor symbol`)) -> df2

write_tsv(df2, 'results/MD2/DiffEdges/upregulatedLR.tsv')

df2 %>% 
  slice_max(`Delta edge specificity weight`, n=50) %>%
  select(`Cell-type pair`, `Ligand-receptor pair`, `Delta edge specificity weight`)

df2 %>% 
  slice_max(`Delta edge specificity weight`, n=50) -> df2_sub

df2_sub %>%
  ggplot(aes(`Cell-type pair`,  `Ligand-receptor pair`)) + 
  geom_raster(aes(fill = `Delta edge specificity weight`)) +
  theme_classic() +
  scale_y_discrete(limits = rev(sort(unique(df2_sub$`Ligand-receptor pair`)))) +
  scale_fill_gradientn(colours=c("lightblue", "darkblue")) +
  theme(
    axis.text.x = element_text(size=5,angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size=5)
    ) +
  labs(caption="only Delta edge specificity weight top 50 pairs")-> g
  
ggsave("results/MD2/DiffEdges/heatmap_upregulated_top50.pdf", g)
