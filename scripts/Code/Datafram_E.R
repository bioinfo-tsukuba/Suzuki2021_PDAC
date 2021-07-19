library(tidyverse)
library(dplyr)

# read Dataframe_D
DFD <- read_csv("results/NATMI_each_patient/Dataframe_AtoE/Dataframe_D.csv")

unique(DFD$cell_type_pair)

# naを0に置換
DFD[is.na(DFD)] <- 0

# make dataframe E
DFD %>%
  group_by(cell_type_pair) %>%
  summarise(mean_NormHRL = mean(NormHRL), mean_NormHRH = mean(NormHRH)) -> DF_E

# save
write_csv(DF_E, "results/NATMI_each_patient/Dataframe_AtoE/Dataframe_E.csv")
