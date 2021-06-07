library(Cairo)
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(circlize)
library(grid)

# read csv of filtered LR pairs, HR<1
df_HMA<- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/NATMI_LR_patients_HMA.csv")

df_HMA %>% rename("Edge_average_expression_weight" = "Edge average expression weight") %>%
  select(LR, cell_type_pair, Edge_average_expression_weight, Patient) %>%
  filter(!is.na(Patient)) ->df_GB

filter(df_GB, is.na(Patient))

# grade 2, P02, P05, P06, P08, P09,P10
subset(df_GB, subset = Patient %in% c( "P02", "P05", "P06", "P08", "P09","P10")) -> df_G2

# without Endo->Endo, Endo->DC, DC->Endo, DC->DC
filter(df_G2, !cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC")) -> df_G2

df_G2 %>%
  group_by(.dots = lapply(c("LR","cell_type_pair"), as.symbol)) %>% 
  summarise(mean = mean(Edge_average_expression_weight)) %>%
  ungroup() ->df_G2_mean

# G2列追加
df_G2_mean <- mutate(df_G2_mean, Grade = "2")

# Check
str(df_G2_mean)

#　重複列がないか調べる
df_G2_mean %>%
  group_by(LR, cell_type_pair) %>%
  filter(n()>1)

# grade 3, P04, P07
subset(df_GB, subset = Patient %in% c( "P04", "P07")) -> df_G3

# without Endo->Endo, Endo->DC, DC->Endo, DC->DC
filter(df_G3, !cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC")) -> df_G3

df_G3 %>%
  group_by(.dots = lapply(c("LR","cell_type_pair"), as.symbol)) %>% 
  summarise(mean = mean(Edge_average_expression_weight)) %>%
  ungroup() ->df_G3_mean

# G2列追加
df_G3_mean <- mutate(df_G3_mean, Grade = "3")

# Check
str(df_G3_mean)

#　重複列がないか調べる
df_G3_mean %>%
  group_by(LR, cell_type_pair) %>%
  filter(n()>1)

# grade 4, P01, P03
subset(df_GB, subset = Patient %in% c( "P01", "P3")) -> df_G4

# without Endo->Endo, Endo->DC, DC->Endo, DC->DC
filter(df_G4, !cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC")) -> df_G4

df_G4 %>%
  group_by(.dots = lapply(c("LR","cell_type_pair"), as.symbol)) %>% 
  summarise(mean = mean(Edge_average_expression_weight)) %>%
  ungroup() ->df_G4_mean

# G2列追加
df_G4_mean <- mutate(df_G4_mean, Grade = "4")

# Check
str(df_G4_mean)

#　重複列がないか調べる
df_G4_mean %>%
  group_by(LR, cell_type_pair) %>%
  filter(n()>1)

# Bind
bind_rows(df_G2_mean, df_G3_mean, df_G4_mean) -> df_result

# rename column
rename(df_result, mean_edge_average_ew = mean) -> df_result

# save df_result
# write_csv(df_result, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/Datafram_B.csv")

# Convert tidy data into matrix
df_result %>% 
  filter(!is.na(Grade)) %>%
  pivot_wider(names_from=Grade, values_from=mean_edge_average_ew, 
              values_fn = mean, values_fill=0) -> df_mat

df_mat %>%
  select(-LR, -cell_type_pair) %>%
  as.matrix() -> mat1

# save df_mat
write_csv(df_mat, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/Datafram_C.csv")

# rownames(mat1) <- df_mat$LR

# Check
str(mat1)

# 図A-1
pdf("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HM_datafram_B.pdf") 

# col = f1をオプションにする
f1 <- colorRamp2(seq(min(mat1), max(mat1), length = 2), c("#EEEEEE", "blue"))

# annotation
ha = rowAnnotation(LR = df_mat$LR, 
                   cell_type_pair = df_mat$cell_type_pair)
                  # labels_gp = gpar(fontsize = 10))

# heatmap
Heatmap(mat1, f1, left_annotation = ha)
       # row_names_gp = gpar(fontsize = 5), column_names_gp =gpar(fontsize = 5))

dev.off() # これでモードが終了し、 PDFができる

# 図A-1 text annotation
pdf("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HM_datafram_B_textanno.pdf",
    width = 10, height = 30) 

# col = f1をオプションにする
f1 <- colorRamp2(seq(min(mat1), max(mat1), length = 2), c("#EEEEEE", "blue"))

# annotation
ha = rowAnnotation(LR = anno_text(df_mat$LR, gp = gpar(fontsize = 2)), 
                   cell_type_pair = anno_text(df_mat$cell_type_pair, gp = gpar(fontsize = 2)))
# labels_gp = gpar(fontsize = 10))

# heatmap
Heatmap(mat1, f1, left_annotation = ha,
        column_names_gp =gpar(fontsize = 5))

dev.off() # これでモードが終了し、 PDFができる


#　図B each cell type

# make directory
dir.create("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HM_Grade/")

ctps <- unique(df_mat$cell_type_pair)
for(ctp in ctps){
  df_mat %>%
    filter(cell_type_pair == ctp) -> df_mat_sub
  
  df_mat_sub %>%
    select(-LR, -cell_type_pair) %>%
    as.matrix() -> mat1_sub
  rownames(mat1_sub) <- df_mat_sub$LR
  pdf_name <- paste0("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HM_Grade/", "heatmap_", gsub("->", "-", ctp), ".pdf")
  pdf(pdf_name)
  ht <- Heatmap(mat1_sub, f1,  name = ctp, 
                row_names_gp = gpar(fontsize = 5), column_names_gp =gpar(fontsize = 5))
  draw(ht)
  dev.off()
}

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
  pdf(pdf_name)
  ht <- Heatmap(mat1_sub, f1, column_order = sort(colnames(mat1_sub)),
                row_order = sort(rownames(mat1_sub)),
                row_names_gp = gpar(fontsize = 5), column_names_gp =gpar(fontsize = 5),
                row_title = ctp )   # Show cell type name on the top
  # Store the heatmap
  list_ht <- append(list_ht, list(ht)) 
}
# PDF
pdf_name <- ("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HMDFB_by_grade.pdf")
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

# 図B TAMだけを抽出する
df_B <- read.csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/Dataframe_AtoE/Datafram_B.csv")

# Convert tidy data into matrix
df_B %>% 
  filter(!is.na(Grade)) %>%
  pivot_wider(names_from=Grade, values_from=mean_edge_average_ew, 
              values_fn = mean, values_fill=0) -> df_mat

df_mat %>%
  select(-LR, -cell_type_pair) %>%
  as.matrix() -> mat1

# Heatmap
ctps <- unique(df_mat$cell_type_pair[str_detect(df_mat$cell_type_pair, "TAM")])

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
pdf_name <- ("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/heatmap_TAM_by_grade.pdf")
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

