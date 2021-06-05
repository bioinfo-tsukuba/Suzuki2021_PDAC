library(Cairo)
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(circlize)
library(grid)

df_mat <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/Dataframe_AtoE/Datafram_C.csv")

dir.create("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HM_Grade_edit/")

ctps <- unique(df_mat$cell_type_pair)
for(ctp in ctps){
  df_mat %>%
    filter(cell_type_pair == ctp) -> df_mat_sub
  
  df_mat_sub %>%
    select(-LR, -cell_type_pair) %>%
    as.matrix() -> mat1_sub
  rownames(mat1_sub) <- df_mat_sub$LR
  pdf_name <- paste0("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HM_Grade_edit/", "heatmap_", gsub("->", "-", ctp), ".pdf")
  pdf(pdf_name)
  f1 <- colorRamp2(seq(min(mat1_sub), max(mat1_sub), length = 2), c("#EEEEEE", "blue"))
  ht <- Heatmap(mat1_sub, f1,  name = ctp, column_order = sort(colnames(mat1_sub)),
                row_order = sort(rownames(mat1_sub)),
                row_names_gp = gpar(fontsize = 5), column_names_gp =gpar(fontsize = 5))
  draw(ht)
  dev.off()
}