###############################################################################
# Setup
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, Seurat)

system("mkdir -p results/MD1")

###############################################################################
# Generate Seurat object
###############################################################################

df_patients <- read_csv("data/WeiLin_pdac10/Patient_data.csv",
  col_names=c("patients","age","gender","diagnosis","Primary_or_Metasitasis","stage","grade"), skip = 1)

df_meta <-
  read_csv( "data/WeiLin_pdac10/PDAC10_meta.csv", col_names=c("cell","patients","type"), skip = 1) %>%
  left_join(df_patients, key = "patients") %>%
  mutate(grade_lowhigh = if_else(grade >2, "high", "low"))

###############################################################################
# Barplot
# x-axis: Patiens, y-axis: % of cell type, facet: cell type
###############################################################################

df_plot <-
  df_meta %>%
  select(patients, type, grade, stage) %>%
  group_by(patients, grade, stage) %>%
  count(type) %>%
  mutate(percent = n/sum(n)*100)

p_grade <-
  ggplot(df_plot, aes(x = patients, y = percent)) +
  geom_col(aes(fill = grade)) +
  ylim(0,100) +
  theme_bw() +
  facet_wrap(~type)

p_stage <-
  ggplot(df_plot, aes(x = patients, y = percent)) +
  geom_col(aes(fill = stage)) +
  ylim(0,100) +
  theme_bw() +
  facet_wrap(~type)

ggsave("results/MD1/percent_celltype_by_grade.pdf", p_grade)
ggsave("results/MD1/percent_celltype_by_stage.pdf", p_stage)


