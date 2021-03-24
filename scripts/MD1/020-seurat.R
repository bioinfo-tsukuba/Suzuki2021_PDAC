###############################################################################
# Setup
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, Seurat, umap, broom, glue)

system("mkdir -p results/MD1")

###############################################################################
# Generate Seurat object
###############################################################################

df_meta <- read_csv( "data/WeiLin_pdac10/PDAC10_meta.csv", col_names=c("cell","id","type"), skip = 1) %>%
  mutate(grade = if_else(id %in% c("P02","P05","P06","P08","P09","P10"), "low", "high"))
head(df_meta)

data_seurat <- CreateSeuratObject(counts = Read10X(data.dir = "data/WeiLin_pdac10/For_Seurat/"), project = "pdac")
data_seurat@meta.data <- cbind(data_seurat@meta.data, df_meta)

###############################################################################
# UMAP by each cell type
###############################################################################

data_umap <- map_dfr(unique(df_meta$type), function(.celltype) {
  .data <- data_seurat[, str_which(df_meta$type, .celltype)]
  .data_umap <-
    FindVariableFeatures(
    object = .data,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

  .data_umap[["umap"]]@"cell.embeddings" %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    as_tibble %>%
    inner_join(df_meta, key = "cell")
})

p_umap <- ggplot(data_umap, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = grade)) +
  theme_bw() +
  facet_wrap(~type, scale = "free")
ggsave("results/MD1/umap_celltype.pdf", p_umap)


