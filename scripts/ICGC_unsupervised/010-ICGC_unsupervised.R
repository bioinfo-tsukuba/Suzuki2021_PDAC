library(Seurat)
library(tidyverse)

paths_icgc_data <- c(
  "data/ICGC/survival_PAAD-US.csv.gz",
  "data/ICGC/survival_PACA-AU.csv.gz",
  "data/ICGC/survival_PACA-CA.csv.gz"
)
outdir = "results/ICGC_unsupervised/20210404"


#####################################
# Make directory
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

#####################################
# Read ICGC data
lapply(paths_icgc_data, function(path_icgc_data){
  read_csv(path_icgc_data) %>%
    mutate(cohort = gsub(".csv.gz", "", gsub("survival_", "", basename(path_icgc_data))))
}) %>%
  bind_rows() -> df1

#########################################
# Make Seurat object
#########################################
# Convert to gene x sample matrix
# [Tentative] Currently, duplicated rows are averaged
df1 %>% 
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
df1 %>%
  select(-exp, -gene) %>%
  distinct() %>%
  right_join(tibble(id=colnames(exp_mat)), by="id") -> df_meta
seurat <- AddMetaData(seurat, metadata=df_meta$sex, col.name = "sex")
seurat <- AddMetaData(seurat, metadata=df_meta$status, col.name = "status")
seurat <- AddMetaData(seurat, metadata=df_meta$time, col.name = "time")
seurat <- AddMetaData(seurat, metadata=df_meta$cohort, col.name = "cohort")

#########################################
# Perform dimensional reduction
#########################################
# Normalize data
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

# Find ~1k variable genes
seurat <- FindVariableFeatures(seurat, selection.method = "vst")

# Scale data
seurat <- ScaleData(seurat)

# Run PCA
seurat <- RunPCA(seurat, verbose = FALSE, approx=FALSE)
ElbowPlot(seurat, ndims = 50)

# Run umap on PCA result
seurat <- RunUMAP(seurat, dims=1:20, seed.use=1234, n.components=2)

#########################################
# Visualize PCA
#########################################
FeaturePlot(seurat, reduction = 'pca', features = "time",
            shape.by = 'sex')
ggsave(file.path(outdir, "pca_by_time_sex.pdf"))

FeaturePlot(seurat, reduction = 'pca', features = "time",
            shape.by = 'cohort')
ggsave(file.path(outdir, "pca_by_time_cohort.pdf"))

DimPlot(seurat, reduction = 'pca', group.by = "cohort",
        shape.by = 'cohort') 
ggsave(file.path(outdir, "pca_by_cohort.pdf"))

#########################################
# Visualize UMAP
#########################################
FeaturePlot(seurat, reduction = 'umap', features = "time",
            shape.by = 'sex') 
ggsave(file.path(outdir, "umap_by_time_sex.pdf"))

FeaturePlot(seurat, reduction = 'umap', features = "time",
            shape.by = 'cohort') 
ggsave(file.path(outdir, "umap_by_time_cohort.pdf"))

DimPlot(seurat, reduction = 'umap', group.by = "cohort",
            shape.by = 'cohort') 
ggsave(file.path(outdir, "umap_by_cohort.pdf"))

#########################################
# Find clusters
#########################################
# Clustering
seurat <- FindNeighbors(seurat, dims = 1:20)
seurat <- FindClusters(seurat, resolution = 0.8)

# Visualize clusters on PCA
DimPlot(seurat, reduction = 'pca', group.by = "seurat_clusters",
        shape.by = 'cohort') 
ggsave(file.path(outdir, "pca_by_cluster_cohort.pdf"))

# Visualize clusters on PCA
DimPlot(seurat, reduction = 'umap', group.by = "seurat_clusters",
        shape.by = 'cohort') 
ggsave(file.path(outdir, "umap_by_cluster_cohort.pdf"))

# Find cluster marker genes
df_seurat_markers <- FindAllMarkers(seurat, min.diff.pct = 0.3, only.pos = TRUE) %>%
  as_tibble()

# Save cluster annotation
tibble(id=names(seurat$seurat_clusters), cluster=seurat$seurat_clusters) -> 
  df_cluster
write_tsv(
  df_cluster, file.path(outdir, "df_cluster.tsv")
)

# Save cluster marker genes
write_tsv(
  df_seurat_markers,
  file.path(outdir, "df_seurat_markers.tsv")
)

#########################################
# Save Seurat object
#########################################
saveRDS(seurat, file = file.path(outdir, "seurat.rds"))

#####################
# sessionInfo()
sessionInfo()
