library(Cairo)
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(circlize)

# read csv
df <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/Fig1/LR_adjPval_meanHR_screened.csv")

# [HM-A] HR<1 の LRペアのヒートマップ
df %>% filter(meanHR <= 1) %>% select(LR) -> HMA

# P01
df_natmi_01 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P01/Edges_lrc2p.csv")
df_01 <- mutate(df_natmi_01, Patient = "P01") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_01, by = "LR") -> df_01
unique(df_01$LR)

# P02
df_natmi_02 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P02/Edges_lrc2p.csv")
df_02 <- mutate(df_natmi_02, Patient = "P02") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_02, by = "LR") -> df_02
unique(df_02$LR)

# P03
df_natmi_03 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P03/Edges_lrc2p.csv")
df_03 <- mutate(df_natmi_03, Patient = "P03") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_03, by = "LR") -> df_03
unique(df_03$LR)

# P04
df_natmi_04 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P04/Edges_lrc2p.csv")
df_04 <- mutate(df_natmi_04, Patient = "P04") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_04, by = "LR") -> df_04
unique(df_04$LR)

# P05
df_natmi_05 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P05/Edges_lrc2p.csv")
df_05 <- mutate(df_natmi_05, Patient = "P05") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_05, by = "LR") -> df_05
unique(df_05$LR)

# P06
df_natmi_06 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P06/Edges_lrc2p.csv")
df_06 <- mutate(df_natmi_06, Patient = "P06") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_06, by = "LR") -> df_06
unique(df_06$LR)

# P07
df_natmi_07 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P07/Edges_lrc2p.csv")
df_07 <- mutate(df_natmi_07, Patient = "P07") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_07, by = "LR") -> df_07
unique(df_07$LR)


# P08
df_natmi_08 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P08/Edges_lrc2p.csv")
df_08 <- mutate(df_natmi_08, Patient = "P08") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_08, by = "LR") -> df_08
unique(df_08$LR)

# P09
df_natmi_09 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P09/Edges_lrc2p.csv")
df_09 <- mutate(df_natmi_09, Patient = "P09") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_09, by = "LR") -> df_09
unique(df_09$LR)

# P10
df_natmi_10 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P10/Edges_lrc2p.csv")
df_10 <- mutate(df_natmi_10, Patient = "P10") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(HMA, df_10, by = "LR") -> df_10
unique(df_10$LR)

# Bind
bind_rows(df_01, df_02, df_03, df_04, df_05, df_06, df_07, df_08, df_09, df_10) -> df_result

write_csv(df_result, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/NATMI_LR_patients_HMA.csv")

# 基本の表
df_result %>% select(LR, cell_type_pair, Patient) -> df_common_LR
unique(df_common_LR$LR)
write_csv(df_common_LR, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/NATMI_common_LR_patients_HMA.csv")

# 集計
df_result %>% 
  filter(!is.na(Patient)) %>%
  group_by(LR, cell_type_pair) %>%
  summarize(number_of_patient = n()) -> df_summary
write_csv(df_summary, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HMA_summary_number_of_patient.csv")
unique(df_summary$LR)

dfdf <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HMA_summary_number_of_patient.csv")
# check how many columns contain "Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC"
filter(dfdf, cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC"))
# without Endo->Endo, Endo->DC, DC->Endo, DC->DC
filter(dfdf, !cell_type_pair %in% c("Endo->Endo", "Endo->DC", "DC->Endo", "DC->DC")) -> df2

# Convert tidy data into matrix
df2 %>% 
  pivot_wider(names_from=cell_type_pair, values_from=number_of_patient, values_fill=0) -> df_mat

df_mat %>%
  select(-LR) %>%
  as.matrix() -> mat1

rownames(mat1) <- df_mat$LR

# Check
str(mat1)

library(circlize)
# この後グラフィックの出力はPDFとして保存しますよモード
pdf("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/HMA_heatmap_all_patients.pdf") 

# col = f1をオプションにする
f1 <- colorRamp2(seq(min(mat1), max(mat1), length = 2), c("#EEEEEE", "blue"))

# heatmap
Heatmap(mat1, col = f1, 
        row_names_gp = gpar(fontsize = 5), column_names_gp =gpar(fontsize = 5))

dev.off() # これでモードが終了し、 PDFができる


