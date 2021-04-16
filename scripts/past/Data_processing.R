#This script is to make the em file for NATMI

library(Seurat)
library(dplyr)
pdac_data <- Read10X(data.dir = "SSD/data/WeiLin_pdac10/For_Seurat/")
pdac <- CreateSeuratObject(counts = pdac_data, project = "pdac")

pdac_expr <- GetAssayData(object = pdac, assay.type = "RNA", slot = "data")
pdac_expr <- as(Class = 'matrix', object = pdac_expr)
write.csv(x = pdac_expr, file = "SSD/processed_data/expression_cluster1.csv", quote = FALSE)