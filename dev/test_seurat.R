library(Seurat)

pdac_data <- Read10X(data.dir = "data/WeiLin_pdac10/For_Seurat/")
pdac <- CreateSeuratObject(counts = pdac_data, project = "pdac")


pdac <- FindVariableFeatures(pdac, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pdac), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pdac)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pdac)
pdac <- ScaleData(pdac, features = all.genes)
pdac <- RunPCA(pdac, features = VariableFeatures(object = pdac))
DimPlot(pdac, reduction = "pca")

df1 <- read.csv("data/WeiLin_pdac10/PDAC10_meta.csv", stringsAsFactors = FALSE)



pdac@meta.data <- cbind(pdac@meta.data, df1[, c("sampleid", "types")])
DimPlot(pdac, reduction = "pca", group.by = "types")
DimPlot(pdac, reduction = "pca", group.by = "sampleid")
