library(tidyverse)

# P01
P01 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P01/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P01 %>% group_by(cell_type_pair) %>% summarise(P01_amount = n_distinct(LR)) -> P01_summary

# P02
P02 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P02/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P02 %>% group_by(cell_type_pair) %>% summarise(P02_amount = n_distinct(LR)) -> P02_summary

# P03
P03 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P03/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P03 %>% group_by(cell_type_pair) %>% summarise(P03_amount = n_distinct(LR)) -> P03_summary

# P04
P04 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P04/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P04 %>% group_by(cell_type_pair) %>% summarise(P04_amount = n_distinct(LR)) -> P04_summary

# P05
P05 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P05/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P05 %>% group_by(cell_type_pair) %>% summarise(P05_amount = n_distinct(LR)) -> P05_summary

# P06
P06 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P06/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P06 %>% group_by(cell_type_pair) %>% summarise(P06_amount = n_distinct(LR)) -> P06_summary

# P07
P07 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P07/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P07 %>% group_by(cell_type_pair) %>% summarise(P07_amount = n_distinct(LR)) -> P07_summary

# P08
P08 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P08/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P08 %>% group_by(cell_type_pair) %>% summarise(P08_amount = n_distinct(LR)) -> P08_summary

# P09
P09 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P09/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P09 %>% group_by(cell_type_pair) %>% summarise(P09_amount = n_distinct(LR)) -> P09_summary

# P10
P10 <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/ExtractEdges/P10/Edges_lrc2p.csv") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->") 

P10 %>% group_by(cell_type_pair) %>% summarise(P10_amount = n_distinct(LR)) -> P10_summary

# join
cell_type_summary <- full_join(P01_summary, P02_summary, by="cell_type_pair") %>% 
  full_join(P03_summary,  by="cell_type_pair") %>%
  full_join(P04_summary,  by="cell_type_pair") %>%
  full_join(P05_summary,  by="cell_type_pair") %>%
  full_join(P06_summary,  by="cell_type_pair") %>%
  full_join(P07_summary,  by="cell_type_pair") %>%
  full_join(P08_summary,  by="cell_type_pair") %>%
  full_join(P09_summary,  by="cell_type_pair") %>%
  full_join(P10_summary,  by="cell_type_pair") 

write.csv(cell_type_summary, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/cell_type_sumamry.csv")
