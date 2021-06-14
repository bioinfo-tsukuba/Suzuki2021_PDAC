library(Cairo)
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(circlize)


df <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/Dataframe_AtoE/Datafram_B.csv")

# Convert tidy data into matrix
df %>% 
  filter(!is.na(Grade)) %>%
  pivot_wider(names_from=Grade, values_from=mean_edge_average_ew, 
              values_fn = mean, values_fill=0) -> df_B

df_B %>% filter(`4` > `2` & `4` > `3`) -> df_mat

write_csv(df_mat, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/Dataframe_AtoE/Datafram_filtered_grade.csv")

df_mat %>%
  select(-LR, -cell_type_pair) %>%
  as.matrix() -> mat1

# 図B-2 cell type pair ordering by patient #
ctps <- unique(df_mat$cell_type_pair)

# List to store all Heatmap object
list_ht <- list()
for(ctp in ctps){
  df_mat %>%
    filter(cell_type_pair == ctp) -> df_mat_sub
  df_mat_sub %>%
    select(-LR, -cell_type_pair) %>%
    as.matrix() -> mat1_sub
  rownames(mat1_sub) <- df_mat_sub$LR
  f1 <- colorRamp2(seq(min(mat1_sub), max(mat1_sub), length = 2), c("#EEEEEE", "blue"))
  ht <- Heatmap(mat1_sub, f1, column_order = sort(colnames(mat1_sub)),
                row_order = sort(rownames(mat1_sub)),
                row_names_gp = gpar(fontsize = 5), column_names_gp =gpar(fontsize = 5),
                row_title = ctp )   # Show cell type name on the top
  # Store the heatmap
  list_ht <- append(list_ht, list(ht)) 
}
# PDF
pdf_name <- ("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/Heatmap_AtoB/HRL_HMgrade_filtered.pdf")
pdf(pdf_name, width = 30, height = 40)
# Determine the numbers of rows and columns
n_plot <- length(list_ht)
n_row <- ceiling(sqrt(n_plot))
n_col <- ceiling(n_plot/n_row)
# Generate n_row x n_col plot subspace
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = n_row, nc = n_col)))
for(i in 1:n_plot){
  # Determine the plot subspace
  current_row <- ceiling(i / n_col)
  current_col <- i - (current_row-1) * n_col
  # Plot heatmap in the subspace
  pushViewport(viewport(layout.pos.row = current_row, layout.pos.col = current_col))
  draw(list_ht[[i]], newpage = FALSE)
  upViewport()
}
upViewport()
dev.off()

