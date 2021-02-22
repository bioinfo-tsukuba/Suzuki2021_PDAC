#heatmap for UP-regulated_mean.csv
library(tidyverse)
path_upregulated <- "results/MD2/DiffEdges/Delta_edges_lrc2p/UP-regulated_mean.csv" #TODO
df1 <- read_csv(path_upregulated, col_names=TRUE)

df1 %>% 
  mutate(`Cell-type pair`=paste0(`Sending cluster`, "->", `Target cluster`)) %>%
  mutate(`Ligand-receptor pair`=paste0(`Ligand symbol`, "->", `Receptor symbol`)) -> df2

write_tsv(df2, 'results/MD2/DiffEdges/tableForHeatmap.tsv')


df2 %>% 
  arrange(-`Delta edge specificity weight`) %>%
  top_n(100) %>%
  ggplot(aes(`Cell-type pair`,  `Ligand-receptor pair`)) + 
  geom_raster(aes(fill = `Delta edge specificity weight`)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=5,angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size=5)
    ) +
  labs(caption="only Delta edge specificity weight top 100 pairs")-> g
  
ggsave("results/MD2/DiffEdges/heatmap_upregulated_top100.pdf",g)


df <- data.frame(
  x = rep(c(2, 5, 7, 9, 12), 2),
  y = rep(c(1, 2), each = 5),
  z = factor(rep(1:5, each = 2)),
  w = rep(diff(c(0, 4, 6, 8, 10, 14)), 2)
)
g <- ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = z), colour = "grey50")
ggsave("g.pdf",g)

