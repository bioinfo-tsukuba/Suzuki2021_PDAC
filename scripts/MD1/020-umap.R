###############################################################################
# Setup
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, Seurat, patchwork)

system("mkdir -p results/MD1")

###############################################################################
# Generate Seurat object
###############################################################################

df_patients <- read_csv("data/WeiLin_pdac10/Patient_data.csv",
  col_names=c("patients","age","gender","diagnosis","Primary_or_Metasitasis","stage","grade"), skip = 1)

df_meta <-
  read_csv( "data/WeiLin_pdac10/PDAC10_meta.csv", col_names=c("cell","patients","type"), skip = 1) %>%
  left_join(df_patients, key = "patients") %>%
  mutate(grade_lowhigh = if_else(grade >2, "high", "low"))

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

p_patients <-
  data_umap %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = patients)) +
  ggtitle("Patients") +
  theme_bw() +
  facet_wrap(~type, scale = "free")

p_grade <-
  data_umap %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = grade)) +
  ggtitle("Grade") +
  theme_bw() +
  facet_wrap(~type, scale = "free")

p_stage <-
  data_umap %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = stage)) +
  ggtitle("Stage") +
  theme_bw() +
  facet_wrap(~type, scale = "free")

ggsave("results/MD1/umap_celltype.pdf",
  (p_patients / p_grade / p_stage),
  width = 13, height = 30)


p_patients_wo2 <-
  data_umap %>%
  filter(grade != 2) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = patients)) +
  ggtitle("patients") +
  theme_bw() +
  facet_wrap(~type, scale = "free")

p_grade_wo2 <-
  data_umap %>%
  filter(grade != 2) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = grade)) +
  ggtitle("grade") +
  theme_bw() +
  facet_wrap(~type, scale = "free")

p_stage_wo2 <-
  data_umap %>%
  filter(grade != 2) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = stage)) +
  ggtitle("Stage") +
  theme_bw() +
  facet_wrap(~type, scale = "free")

ggsave("results/MD1/umap_celltype_wo_g2.pdf",
  (p_patients_wo2 / p_grade_wo2 / p_stage_wo2),
  width = 13, height = 30)
