library(tidyverse)
library(dplyr)

# HR<1 = HRL
HRL <- read_csv("results/NATMI_each_patient/NATMI_LR_patients_HMA.csv")
unique(HRL$LR)

# HR>1 = HRH
HRH <- read_csv("results/NATMI_each_patient/NATMI_LR_patients_HMB.csv")
unique(HRH$LR)

# P01 
df_natmi_01 <- read_csv("results/NATMI_each_patient/ExtractEdges/P01/Edges_lrc2p.csv")
df_01 <- mutate(df_natmi_01, Patient = "P01") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_01 %>% select(cell_type_pair, LR, Patient) -> P01
P01 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P01

str(P01)

# HRL
HRL[HRL$Patient=="P01",] -> HRL_P01
unique(HRL_P01$LR)

HRL_P01 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P01

# HRH
HRH[HRH$Patient=="P01",] -> HRH_P01
unique(HRH_P01$LR)

HRH_P01 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P01

# Join P01
full_join(P01, HRL_P01, by = "cell_type_pair") -> P01
full_join(P01, HRH_P01, by = "cell_type_pair") -> P01

mutate(P01, Patient = "P01") -> P01

# P02
df_natmi_02 <- read_csv("results/NATMI_each_patient/ExtractEdges/P02/Edges_lrc2p.csv")
df_02 <- mutate(df_natmi_02, Patient = "P02") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_02 %>% select(cell_type_pair, LR, Patient) -> P02
P02 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P02

str(P02)

# HRL
HRL[HRL$Patient=="P02",] -> HRL_P02
unique(HRL_P02$LR)

HRL_P02 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P02

# HRH
HRH[HRH$Patient=="P02",] -> HRH_P02
unique(HRH_P02$LR)

HRH_P02 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P02

# Join P02
full_join(P02, HRL_P02, by = "cell_type_pair") -> P02
full_join(P02, HRH_P02, by = "cell_type_pair") -> P02

mutate(P02, Patient = "P02") -> P02

# P03
df_natmi_03 <- read_csv("results/NATMI_each_patient/ExtractEdges/P03/Edges_lrc2p.csv")
df_03 <- mutate(df_natmi_03, Patient = "P03") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_03 %>% select(cell_type_pair, LR, Patient) -> P03
P03 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P03

str(P03)

# HRL
HRL[HRL$Patient=="P03",] -> HRL_P03
unique(HRL_P03$LR)

HRL_P03 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P03

# HRH
HRH[HRH$Patient=="P03",] -> HRH_P03
unique(HRH_P03$LR)

HRH_P03 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P03

# Join P03
full_join(P03, HRL_P03, by = "cell_type_pair") -> P03
full_join(P03, HRH_P03, by = "cell_type_pair") -> P03

mutate(P03, Patient = "P03") -> P03

# P04
df_natmi_04 <- read_csv("results/NATMI_each_patient/ExtractEdges/P04/Edges_lrc2p.csv")
df_04 <- mutate(df_natmi_04, Patient = "P04") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_04 %>% select(cell_type_pair, LR, Patient) -> P04
P04 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P04

str(P04)

# HRL
HRL[HRL$Patient=="P04",] -> HRL_P04
unique(HRL_P04$LR)

HRL_P04 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P04

# HRH
HRH[HRH$Patient=="P04",] -> HRH_P04
unique(HRH_P04$LR)

HRH_P04 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P04

# Join P04
full_join(P04, HRL_P04, by = "cell_type_pair") -> P04
full_join(P04, HRH_P04, by = "cell_type_pair") -> P04

mutate(P04, Patient = "P04") -> P04

# P05
df_natmi_05 <- read_csv("results/NATMI_each_patient/ExtractEdges/P05/Edges_lrc2p.csv")
df_05 <- mutate(df_natmi_05, Patient = "P05") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_05 %>% select(cell_type_pair, LR, Patient) -> P05
P05 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P05

str(P05)

# HRL
HRL[HRL$Patient=="P05",] -> HRL_P05
unique(HRL_P05$LR)

HRL_P05 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P05

# HRH
HRH[HRH$Patient=="P05",] -> HRH_P05
unique(HRH_P05$LR)

HRH_P05 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P05

# Join P05
full_join(P05, HRL_P05, by = "cell_type_pair") -> P05
full_join(P05, HRH_P05, by = "cell_type_pair") -> P05

mutate(P05, Patient = "P05") -> P05

# P06
df_natmi_06 <- read_csv("results/NATMI_each_patient/ExtractEdges/P06/Edges_lrc2p.csv")
df_06 <- mutate(df_natmi_06, Patient = "P06") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_06 %>% select(cell_type_pair, LR, Patient) -> P06
P06 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P06

str(P06)

# HRL
HRL[HRL$Patient=="P06",] -> HRL_P06
unique(HRL_P06$LR)

HRL_P06 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P06

# HRH
HRH[HRH$Patient=="P06",] -> HRH_P06
unique(HRH_P06$LR)

HRH_P06 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P06

# Join P06
full_join(P06, HRL_P06, by = "cell_type_pair") -> P06
full_join(P06, HRH_P06, by = "cell_type_pair") -> P06

mutate(P06, Patient = "P06") -> P06

# P07
df_natmi_07 <- read_csv("results/NATMI_each_patient/ExtractEdges/P07/Edges_lrc2p.csv")
df_07 <- mutate(df_natmi_07, Patient = "P07") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_07 %>% select(cell_type_pair, LR, Patient) -> P07
P07 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P07

str(P07)

# HRL
HRL[HRL$Patient=="P07",] -> HRL_P07
unique(HRL_P07$LR)

HRL_P07 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P07

# HRH
HRH[HRH$Patient=="P07",] -> HRH_P07
unique(HRH_P07$LR)

HRH_P07 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P07

# Join P07
full_join(P07, HRL_P07, by = "cell_type_pair") -> P07
full_join(P07, HRH_P07, by = "cell_type_pair") -> P07

mutate(P07, Patient = "P07") -> P07

# P08
df_natmi_08 <- read_csv("results/NATMI_each_patient/ExtractEdges/P08/Edges_lrc2p.csv")
df_08 <- mutate(df_natmi_08, Patient = "P08") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_08 %>% select(cell_type_pair, LR, Patient) -> P08
P08 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P08

str(P08)

# HRL
HRL[HRL$Patient=="P08",] -> HRL_P08
unique(HRL_P08$LR)

HRL_P08 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P08

# HRH
HRH[HRH$Patient=="P08",] -> HRH_P08
unique(HRH_P08$LR)

HRH_P08 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P08

# Join P08
full_join(P08, HRL_P08, by = "cell_type_pair") -> P08
full_join(P08, HRH_P08, by = "cell_type_pair") -> P08

mutate(P08, Patient = "P08") -> P08

# P09
df_natmi_09 <- read_csv("results/NATMI_each_patient/ExtractEdges/P09/Edges_lrc2p.csv")
df_09 <- mutate(df_natmi_09, Patient = "P09") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_09 %>% select(cell_type_pair, LR, Patient) -> P09
P09 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P09

str(P09)

# HRL
HRL[HRL$Patient=="P09",] -> HRL_P09
unique(HRL_P09$cell_type_pair)

HRL_P09 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P09

# HRH
HRH[HRH$Patient=="P09",] -> HRH_P09
unique(HRH_P09$cell_type_pair)

HRH_P09 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P09

# Join P09
full_join(P09, HRL_P09, by = "cell_type_pair") -> P09
full_join(P09, HRH_P09, by = "cell_type_pair") -> P09

mutate(P09, Patient = "P09") -> P09

# P10
df_natmi_10 <- read_csv("results/NATMI_each_patient/ExtractEdges/P10/Edges_lrc2p.csv")
df_10 <- mutate(df_natmi_10, Patient = "P10") %>%
  unite("LR", c("Ligand symbol", "Receptor symbol"), sep = "->") %>%
  unite("cell_type_pair", c("Sending cluster", "Target cluster"), sep = "->")

df_10 %>% select(cell_type_pair, LR, Patient) -> P10
P10 %>% 
  group_by(cell_type_pair) %>%
  summarize(NATMI_LR_ALL = n()) -> P10

str(P10)

# HRL
HRL[HRL$Patient=="P10",] -> HRL_P10
unique(HRL_P10$cell_type_pair)

HRL_P10 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRL = n()) -> HRL_P10

# HRH
HRH[HRH$Patient=="P10",] -> HRH_P10
unique(HRH_P10$cell_type_pair)

HRH_P10 %>%
  filter(!is.na(Patient)) %>%
  select(cell_type_pair, LR, Patient) %>%
  group_by(cell_type_pair) %>%
  summarize(HRH = n()) -> HRH_P10

# Join P10
full_join(P10, HRL_P10, by = "cell_type_pair") -> P10
full_join(P10, HRH_P10, by = "cell_type_pair") -> P10

mutate(P10, Patient = "P10") -> P10

# Bind
bind_rows(P01, P02, P03, P04, P05, P06, P07, P08, P09, P10) -> DFD

# Replace NA in HRL and HRH columns to 0
DFD %>%
  replace_na(list(HRL=0, HRH=0)) -> DFD

# row HRL/all
mutate(DFD, NormHRL = (HRL)/(NATMI_LR_ALL)) -> DFD

# row HRH/all
mutate(DFD, NormHRH = (HRH)/(NATMI_LR_ALL)) -> DFD

# save
write_csv(DFD, "results/NATMI_each_patient/Dataframe_AtoE/Dataframe_D.csv")
