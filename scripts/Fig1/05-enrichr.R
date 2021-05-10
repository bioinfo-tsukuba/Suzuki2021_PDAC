################################################################################
# Initialization
################################################################################

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, enrichR)

################################################################################
# Input and format
################################################################################

df <- read_csv("results/Fig1/LR_selected.csv") %>%
  mutate(lr_pair = LR) %>%
  separate(lr_pair, c("ligand", "receptor"), sep = "->") %>%
  pivot_longer(c(ligand, receptor), names_to = "lr", values_to = "gene") %>%
  select(lr, gene, HR) %>%
  distinct()


dbs <- listEnrichrDbs()

# Gene Ontology 2018
dbs_go <- dbs %>%
  select(libraryName) %>%
  filter(str_detect(libraryName, "GO")) %>%
  tail(3) %>%
  pull()

# Reactome, WikiPathways, KEGG
dbs_path <- dbs %>%
  select(libraryName) %>%
  filter(str_detect(libraryName, "Reactome|_Human$")) %>%
  tail(3) %>%
  pull()

dbs_cat <- append(dbs_go, dbs_path)

df_enriched <- unique(df$HR) %>%
  map_dfr(function(.x) {
    map_dfr(dbs_cat, function(.dbs) {
      df %>%
        filter(HR == .x) %>%
        pull(gene) %>%
        enrichr(.dbs) %>%
        flatten_dfr() %>%
        mutate(DB = .dbs) %>%
        mutate(HR = .x)
    })
  }) %>%
  mutate(minusLogAdjPval = -log10(Adjusted.P.value)) %>%
  select(HR, DB, Term, minusLogAdjPval, Genes)

p_col <- df_enriched %>%
  group_by(DB) %>%
  slice_max(minusLogAdjPval, n = 10) %>%
  mutate(Term = reorder(Term, minusLogAdjPval)) %>%
  ggplot(aes(x = minusLogAdjPval, y = Term)) +
  geom_col() +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  facet_wrap(~ HR + DB, scale = "free", ncol = 3)

ggsave("results/Fig1/enrichments.pdf", p_col, width = 30, height = 10)
