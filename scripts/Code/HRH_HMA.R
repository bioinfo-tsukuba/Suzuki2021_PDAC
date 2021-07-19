library(Cairo)
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(circlize)
library(grid)

# read csv of filtered LR pairs, HR>1
df_HMA <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HRH/NATMI_LR_patients_HRH.csv")

df_HMA %>% rename("Edge_average_expression_weight" = "Edge average expression weight") -> df_HMA


# 基本の表+weight
df_HMA %>% select(LR, cell_type_pair, Edge_average_expression_weight, Patient) -> df_common_LR
unique(df_common_LR$LR)
write_csv(df_common_LR, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HRH/HRH_with_weight.csv")

filter(df_common_LR, cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC"))
# without Endo->Endo, Endo->DC, DC->Endo, DC->DC
filter(df_common_LR, !cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC")) -> df1

# Convert tidy data into matrix
df1 %>%
  filter(!is.na(Patient)) %>%
  pivot_wider(
    names_from = Patient, values_from = Edge_average_expression_weight,
    values_fn = mean, values_fill = 0
  ) -> df_mat

df_mat %>%
  select(-LR, -cell_type_pair) %>%
  as.matrix() -> mat1

rownames(mat1) <- df_mat$LR

# Check
str(mat1)

# 図A-1
pdf("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HRH/HRH_HM/HRH_HMA.pdf")

# col = f1をオプションにする
f1 <- colorRamp2(seq(min(mat1), max(mat1), length = 2), c("#EEEEEE", "blue"))

# annotation
ha <- rowAnnotation(
  LR = df_mat$LR,
  cell_type_pair = df_mat$cell_type_pair
)
# labels_gp = gpar(fontsize = 10))

# heatmap
Heatmap(mat1, f1,
  left_annotation = ha,
  row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5)
)

dev.off() # これでモードが終了し、 PDFができる

# 図A-1 text annotation
pdf("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HRH/HRH_HM/HRH_HMA_textanno.pdf",
  width = 10, height = 30
)

# col = f1をオプションにする
f1 <- colorRamp2(seq(min(mat1), max(mat1), length = 2), c("#EEEEEE", "blue"))

# annotation
ha <- rowAnnotation(
  LR = anno_text(df_mat$LR, gp = gpar(fontsize = 2)),
  cell_type_pair = anno_text(df_mat$cell_type_pair, gp = gpar(fontsize = 2))
)
# labels_gp = gpar(fontsize = 10))

# heatmap
Heatmap(mat1, f1,
  left_annotation = ha,
  column_names_gp = gpar(fontsize = 5)
)

dev.off() # これでモードが終了し、 PDFができる


# split the heatmap 図A-2
pdf("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HRH/HRH_HM/HRH_HMA_split.pdf")

# col = f1をオプションにする
f1 <- colorRamp2(seq(min(mat1), max(mat1), length = 2), c("#EEEEEE", "blue"))

# heatmap
Heatmap(mat1, f1,
  row_split = df_mat$cell_type_pair,
  row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5)
)

dev.off() # これでモードが終了し、 PDFができる

# make directory
dir.create("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HRH/sub_set/")

# A-2
ctps <- unique(df_mat$cell_type_pair)
for (ctp in ctps) {
  df_mat %>%
    filter(cell_type_pair == ctp) -> df_mat_sub

  df_mat_sub %>%
    select(-LR, -cell_type_pair) %>%
    as.matrix() -> mat1_sub
  rownames(mat1_sub) <- df_mat_sub$LR
  pdf_name <- paste0("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HRH/sub_set/", "heatmap_", gsub("->", "-", ctp), ".pdf")
  pdf(pdf_name)
  f1 <- colorRamp2(seq(min(mat1_sub), max(mat1_sub), length = 2), c("#EEEEEE", "blue"))
  ht <- Heatmap(mat1_sub, f1,
    name = ctp,
    column_order = sort(colnames(mat1_sub)),
    row_order = sort(rownames(mat1_sub)),
    row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5)
  )
  draw(ht)
  dev.off()
}

# 図A-2 df_1のcell type pair ordering by patient #
ctps <- unique(df_mat$cell_type_pair)

# List to store all Heatmap object
list_ht <- list()
for (ctp in ctps) {
  df_mat %>%
    filter(cell_type_pair == ctp) -> df_mat_sub
  df_mat_sub %>%
    select(-LR, -cell_type_pair) %>%
    as.matrix() -> mat1_sub
  rownames(mat1_sub) <- df_mat_sub$LR
  f1 <- colorRamp2(seq(min(mat1_sub), max(mat1_sub), length = 2), c("#EEEEEE", "blue"))
  ht <- Heatmap(mat1_sub, f1,
    column_order = sort(colnames(mat1_sub)),
    row_order = sort(rownames(mat1_sub)),
    row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
    row_title = ctp
  ) # Show cell type name on the top
  # Store the heatmap
  list_ht <- append(list_ht, list(ht))
}
# PDF
pdf_name <- ("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HRH/HRH_HM_each_cellTypePair_ordered.pdf")
pdf(pdf_name, width = 30, height = 40)
# Determine the numbers of rows and columns
n_plot <- length(list_ht)
n_row <- ceiling(sqrt(n_plot))
n_col <- ceiling(n_plot / n_row)
# Generate n_row x n_col plot subspace
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = n_row, nc = n_col)))
for (i in 1:n_plot) {
  # Determine the plot subspace
  current_row <- ceiling(i / n_col)
  current_col <- i - (current_row - 1) * n_col
  # Plot heatmap in the subspace
  pushViewport(viewport(layout.pos.row = current_row, layout.pos.col = current_col))
  draw(list_ht[[i]], newpage = FALSE)
  upViewport()
}
upViewport()
dev.off()
