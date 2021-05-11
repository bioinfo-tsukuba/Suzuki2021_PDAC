# Load library
library(Seurat)
library(tidyverse)

df_survival <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/survival_meta_analysis/NATMI_meta_filtered_Storey.csv")
df_survival_LR <- select(.data = df_survival, LR)
df_natmi_01 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P01/Edges_lrc2p.csv")

# P01
df_01 <- mutate(df_natmi_01, Patient = "P01") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_01, by = "LR") -> df_01

# P02
df_natmi_02 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P02/Edges_lrc2p.csv")
df_02 <- mutate(df_natmi_02, Patient = "P02") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_02, by = "LR") -> df_02

# P03
df_natmi_03 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P03/Edges_lrc2p.csv")
df_03 <- mutate(df_natmi_03, Patient = "P03") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_03, by = "LR") -> df_03


# P04
df_natmi_04 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P04/Edges_lrc2p.csv")
df_04 <- mutate(df_natmi_04, Patient = "P04") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_04, by = "LR") -> df_04

# P05
df_natmi_05 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P05/Edges_lrc2p.csv")
df_05 <- mutate(df_natmi_05, Patient = "P05") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_05, by = "LR") -> df_05

# P06
df_natmi_06 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P06/Edges_lrc2p.csv")
df_06 <- mutate(df_natmi_06, Patient = "P06") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_06, by = "LR") -> df_06

# P07
df_natmi_07 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P07/Edges_lrc2p.csv")
df_07 <- mutate(df_natmi_07, Patient = "P07") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_07, by = "LR") -> df_07

# P08
df_natmi_08 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P08/Edges_lrc2p.csv")
df_08 <- mutate(df_natmi_08, Patient = "P08") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_08, by = "LR") -> df_08

# P09
df_natmi_09 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P09/Edges_lrc2p.csv")
df_09 <- mutate(df_natmi_09, Patient = "P09") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_09, by = "LR") -> df_09

# P10
df_natmi_10 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P10/Edges_lrc2p.csv")
df_10 <- mutate(df_natmi_10, Patient = "P10") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

left_join(df_survival_LR, df_10, by = "LR") -> df_10

# Bind
bind_rows(df_01, df_02, df_03, df_04, df_05, df_06, df_07, df_08, df_09, df_10) -> df_result

write_csv(df_result, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/NATMI_LR_patients_sumamry.csv")

# 基本の表
select(.data = df_result, LR, cell_type_pair, Patient) -> df_natmi_common_LR_patients
write_csv(df_natmi_common_LR_patients, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/NATMI_common_LR_patients.csv")
