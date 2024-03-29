# This script makes NATMI input (emFile and annFile) for each sample (patient)

# Definition of input and output paths
path_10xdata <- "data/WeiLin_pdac10/For_Seurat/"
path_cell_metadata <- "data/WeiLin_pdac10/PDAC10_meta.csv"
path_outdir <- "processed_data/MD2"

# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Load library
library(Seurat)
library(tidyverse)


# Define P grouping
P_low_grade <- c(
  "P02",
  "P05",
  "P06",
  "P08",
  "P09",
  "P10"
)

P_high_grade <- c(
  "P01",
  "P03",
  "P04",
  "P07"
)

list_P_grouping <- list(P_low_grade = P_low_grade, P_high_grade = P_high_grade)

# Read data
pdac_data <- Read10X(data.dir = path_10xdata)

# Create a Seurat object
pdac <- CreateSeuratObject(counts = pdac_data, project = "pdac")

# Read cell metadata
df_meta <- read_csv(path_cell_metadata, col_names = TRUE)

# Check input
head(df_meta)

# Add cell metadata to the Seurat object
pdac@meta.data <- cbind(pdac@meta.data, df_meta[, c("sampleid", "types")])

# Define sample_id
sample_ids <- unique(df_meta$sampleid)

# For each sampleid
for (i in 1:length(list_P_grouping)) {
  id <- names(list_P_grouping)[i]
  patients <- list_P_grouping[[i]]


  # Get subset of cell names
  df_meta %>%
    filter(sampleid %in% patients) %>%
    .$X1 -> cell_subset

  # Get indicaters of cells from a specific sample (patient)
  cells_subset_seurat <- WhichCells(object = pdac, cells = cell_subset)

  # Extract gene expression data of cells from the specific sample (patient)
  pdac_expr_subset <- GetAssayData(object = pdac, assay = "RNA", slot = "data")[, cells_subset_seurat]

  # Convert subsetted assay data to matrix
  pdac_expr_subset <- as(Class = "matrix", object = pdac_expr_subset)

  # Save expression matrixa as csv
  path_output_emFile <- file.path(path_outdir, paste0(id, "_expression.csv"))
  write.csv(x = pdac_expr_subset, file = path_output_emFile, quote = FALSE)

  # Make metadata for cells from the specific sample (patient)
  pdac_meta_subset <- pdac@meta.data[cells_subset_seurat, ] %>%
    as_tibble(rownames = "barcode") %>%
    rename(annotation = types) %>%
    select(barcode, annotation)

  # Save metada as csv
  path_output_annFile <- file.path(path_outdir, paste0(id, "_annot.csv"))
  write_csv(pdac_meta_subset, path_output_annFile)
}
