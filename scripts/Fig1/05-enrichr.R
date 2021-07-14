################################################################################
# Initialization
################################################################################

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tidytext, enrichR)

################################################################################
# Input and format
################################################################################

df <- read_csv("results/Fig1/LR_adjPval_meanHR_screened.csv") %>%
  mutate(HR = if_else(meanHR > 1, "HR>1", "HR<1")) %>%
  mutate(lr_pair = LR) %>%
  separate(lr_pair, c("ligand", "receptor"), sep = "->") %>%
  pivot_longer(c(ligand, receptor), names_to = "lr", values_to = "gene") %>%
  select(lr, gene, HR) %>%
  distinct()

dbs_all <- listEnrichrDbs()

# Gene Ontology 2018
dbs_go <- dbs_all %>%
  select(libraryName) %>%
  filter(str_detect(libraryName, "GO")) %>%
  tail(3) %>%
  pull()

# BioPlanet, Reactome, WikiPathways, KEGG
dbs_path <- dbs_all %>%
  select(libraryName) %>%
  filter(str_detect(libraryName, "BioPlanet|Reactome|_Human$")) %>%
  tail(4) %>%
  pull()

dbs <- append(dbs_go, dbs_path)

df_enriched <-
  map_dfr(unique(df$HR), function(.hr) {
    map_dfr(unique(df$lr), function(.lr) {
      map_dfr(dbs, function(.dbs) {
        df %>%
          filter(lr == .lr) %>%
          filter(HR == .hr) %>%
          pull(gene) %>%
          enrichr(.dbs) %>%
          flatten_dfr() %>%
          mutate(DB = .dbs) %>%
          mutate(lr = .lr) %>%
          mutate(HR = .hr)
      })
    })
  }) %>%
  mutate(minusLogAdjPval = -log10(Adjusted.P.value)) %>%
  select(lr, HR, DB, Term, minusLogAdjPval, Genes)

p_col <- df_enriched %>%
  group_by(lr, HR, DB) %>%
  slice_max(minusLogAdjPval, n = 10) %>%
  mutate(Term = reorder_within(Term, minusLogAdjPval, lr)) %>%
  mutate(Term = reorder_within(Term, minusLogAdjPval, HR)) %>%
  mutate(Term = reorder_within(Term, minusLogAdjPval, DB)) %>%
  ggplot(aes(x = minusLogAdjPval, y = Term, fill = "FF99CC")) +
  geom_col() +
  theme_bw() +
  scale_x_reverse() +
  theme(text = element_text(size = 15)) +
  facet_wrap(vars(HR, lr, DB), scale = "free", ncol = 4)

p_tile <- df_enriched %>%
  group_by(lr, HR, DB) %>%
  slice_max(minusLogAdjPval, n = 10) %>%
  mutate(Term = reorder_within(Term, minusLogAdjPval, lr)) %>%
  mutate(Term = reorder_within(Term, minusLogAdjPval, HR)) %>%
  mutate(Term = reorder_within(Term, minusLogAdjPval, DB)) %>%
  separate_rows(Genes, sep = ";", convert = TRUE) %>%
  ggplot(aes(x = Genes, y = Term, fill = "FF99CC")) +
  geom_tile() +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  facet_wrap(vars(HR, lr, DB), scale = "free", ncol = 4)


ggsave("results/Fig1/enrichr_barplot.pdf", p_col, width = 70, height = 40, limitsize = FALSE)
ggsave("results/Fig1/enrichr_tile.pdf", p_tile, width = 70, height = 40, limitsize = FALSE)
