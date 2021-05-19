# remove.packages("Cairo")
# install.packages("Cairo")
library(Cairo)
library(ComplexHeatmap)
library(tidyverse)

# read csv
df <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/summary_number_of_patient.csv")


# Convert tidy data into matrix
df %>% 
  pivot_wider(names_from=cell_type_pair, values_from=number_of_patient, values_fill=0) -> df_mat1
df_mat1 %>%
  select(-LR) %>%
  as.matrix() -> mat1
rownames(mat1) <- df_mat1$LR
# Check
str(mat1)

library(circlize)
# この後グラフィックの出力はPDFとして保存しますよモード
pdf("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/heatmap_all_patients.pdf") 

# col = f1をオプションにする
f1 <- colorRamp2(seq(min(mat1), max(mat1), length = 2), c("#EEEEEE", "blue"))

# heatmap
Heatmap(mat1, col = f1, 
        row_names_gp = gpar(fontsize = 5), column_names_gp =gpar(fontsize = 5))

dev.off() # これでモードが終了し、 PDFができる

# check how many columns contain "Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC"
filter(df, cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC"))
# without Endo->Endo, Endo->DC, DC->Endo, DC->DC
filter(df, !cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC")) -> df2

# Convert tidy data into matrix
df2 %>% 
  pivot_wider(names_from=cell_type_pair, values_from=number_of_patient, values_fill=0) -> df_mat2
df_mat2 %>%
  select(-LR) %>%
  as.matrix() -> mat2
rownames(mat2) <- df_mat2$LR
# Check
str(mat2)

#heatmap
pdf("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/heatmap_all_patients_withoutENDO_DC.pdf") 

# col = f1をオプションにする
f1 <- colorRamp2(seq(min(mat2), max(mat2), length = 2), c("#EEEEEE", "blue"))

# Complex heatmap
Heatmap(mat2, col = f1, 
        row_names_gp = gpar(fontsize = 5), column_names_gp =gpar(fontsize = 5))

dev.off() # これでモードが終了し、 PDFができる
