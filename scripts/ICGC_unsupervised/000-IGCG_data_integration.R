# [Note] This script follows a tutorial at https://satijalab.org/seurat/articles/integration_introduction.html.

library(Seurat)
library(tidyverse)
library(patchwork)

paths_icgc_data <- c(
  "data/ICGC/survival_PAAD-US.csv.gz",
  "data/ICGC/survival_PACA-AU.csv.gz",
  "data/ICGC/survival_PACA-CA.csv.gz"
)
outdir <- "results/ICGC_unsupervised/000-IGCG_data_integration"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  outdir <- args[1]
}

#####################################
# Make directory
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


#####################################
# Define function
# Convert to gene x sample matrix
# [Tentative] Currently, duplicated rows are averaged
makeSeuratObject <- function(df_long) {
  df_long %>%
    group_by(gene, id) %>%
    summarise(exp = mean(exp)) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = gene,
      names_from = id,
      values_from = exp
    ) -> df_tmp

  df_tmp %>%
    select(-gene) %>%
    as.matrix() -> exp_mat

  rownames(exp_mat) <- df_tmp$gene

  ###############
  # [Tentative] The reason is not clear but there were many NA values. These values were converted to zero.
  exp_mat[is.na(exp_mat)] <- 0

  # Make Seurat object
  seurat <- CreateSeuratObject(counts = exp_mat)

  # Add metadata
  df_long %>%
    select(-exp, -gene) %>%
    distinct() %>%
    left_join(tibble(id = colnames(exp_mat), num = 1:ncol(exp_mat)), by = "id") %>%
    arrange(num) %>%
    select(-num) -> df_meta

  if (!all(df_meta$id == colnames(exp_mat))) {
    simpleError("IDs are not correctly ordered")
  }

  seurat <- AddMetaData(seurat, metadata = df_meta$sex, col.name = "sex")
  seurat <- AddMetaData(seurat, metadata = df_meta$status, col.name = "status")
  seurat <- AddMetaData(seurat, metadata = df_meta$time, col.name = "time")
  seurat <- AddMetaData(seurat, metadata = df_meta$cohort, col.name = "cohort")

  return(seurat)
}



#####################################
# Read ICGC data and create seurat objects
lapply(paths_icgc_data, function(path_icgc_data) {
  read_csv(path_icgc_data) %>%
    mutate(cohort = gsub(".csv.gz", "", gsub("survival_", "", basename(path_icgc_data)))) %>%
    makeSeuratObject()
}) -> list_seurat_raw

#########################################
# normalize and identify variable features for each dataset independently
list_seurat <- lapply(X = list_seurat_raw, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst")
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list_seurat)

# identify anchors for the use of integrating the datasets together with IntegrateData().
# [Note] Lowered `k.filter`
anchors <- FindIntegrationAnchors(object.list = list_seurat, anchor.features = features, k.filter = 30)

# Creates an 'integrated' data assay
# [Note] Lowered k.weight to avoid an error. See https://github.com/satijalab/seurat/issues/3930.
seurat_combined <- IntegrateData(anchorset = anchors, k.weight = 30)

# specify that we will perform downstream analysis on the corrected data note that the original
DefaultAssay(seurat_combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat_combined <- ScaleData(seurat_combined, verbose = FALSE)
seurat_combined <- RunPCA(seurat_combined, verbose = FALSE)
seurat_combined <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:20)
seurat_combined <- FindNeighbors(seurat_combined, reduction = "pca", dims = 1:20)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(seurat_combined, reduction = "umap", group.by = "cohort")
p2 <- DimPlot(seurat_combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
ggsave(file.path(outdir, "umap_cohort_cluster.svg"), width = 12, height = 9)


p00 <- DimPlot(seurat_combined, reduction = "pca", group.by = "cohort")
p01 <- DimPlot(seurat_combined, reduction = "pca", label = TRUE, repel = TRUE)
p00 + p01
ggsave(file.path(outdir, "pca_cohort_cluster.svg"), width = 12, height = 9)

ElbowPlot(seurat_combined, ndims = 50) + theme(aspect.ratio = 1)
ggsave(file.path(outdir, "pca_elbow.svg"), width = 12, height = 12)

# Save cluster annotation
tibble(
  id = names(seurat_combined$seurat_clusters),
  cluster = seurat_combined$seurat_clusters,
  cohort = seurat_combined$cohort,
  sex = seurat_combined$sex,
  status = seurat_combined$status,
  time = seurat_combined$time,
) -> df_cluster
write_tsv(
  df_cluster, file.path(outdir, "df_cluster.tsv")
)

#####################################
# Find markers
#####################################
# Find cluster marker genes
df_seurat_markers <- FindAllMarkers(seurat_combined, assay = "RNA", min.diff.pct = 0.1, only.pos = TRUE) %>%
  as_tibble()

FeaturePlot(seurat_combined,
  features = df_seurat_markers$gene,
  reduction = "umap", split.by = "seurat_clusters", max.cutoff = 3,
  cols = c("grey", "red")
)
ggsave(file.path(outdir, "umap_marker_genes.svg"), width = 4, height = 24)

# Save cluster marker genes
write_tsv(
  df_seurat_markers,
  file.path(outdir, "df_seurat_markers.tsv")
)

# Violin plot
VlnPlot(seurat_combined,
  features = df_seurat_markers$gene,
  group.by = "seurat_clusters"
) + theme(aspect.ratio = 1)
ggsave(file.path(outdir, "violin_markers.svg"), width = 12, height = 12)

#########################################
#
path_survival_analysis <- "results/survival_meta_analysis/20210417-115137-ab27039/results_LR_meta.csv"
df_survival <- read_csv(path_survival_analysis)

# Filter survival data
## q-value (Storey) < 0.1 && HR < 1 in all cohorts
df_survival %>%
  separate(col = LR, into = c("ligand", "receptor"), remove = FALSE, sep = "->") %>%
  filter(meta_qval_storey < 0.1, `hr_PAAD-US` < 1, `hr_PACA-AU` < 1, `hr_PACA-CA` < 1) -> df_survival

for (i in 1:nrow(df_survival)) {
  l <- df_survival[i, ]$ligand
  r <- df_survival[i, ]$receptor
  FeaturePlot(seurat_combined,
    reduction = "umap", features = c(l, r),
    shape.by = "cohort"
  ) + theme(aspect.ratio = 1)
  ggsave(file.path(outdir, paste0("umap_", l, "_", r, ".pdf")))
}

for (i in 1:nrow(df_survival)) {
  l <- df_survival[i, ]$ligand
  r <- df_survival[i, ]$receptor
  VlnPlot(seurat_combined,
    features = c(l, r),
    group.by = "seurat_clusters"
  ) + theme(aspect.ratio = 1)
  ggsave(file.path(outdir, paste0("vln_", l, "_", r, ".pdf")))
}


#########################################
# Save Seurat object
#########################################
saveRDS(seurat_combined, file = file.path(outdir, "seurat_combined.rds"))

#####################################
# sessionInfo()
sessionInfo()
