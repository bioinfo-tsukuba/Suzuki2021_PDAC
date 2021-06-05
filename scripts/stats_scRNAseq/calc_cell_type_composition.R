path_data <- "data/WeiLin_pdac10/PDAC10_meta.csv"
path_outdir <- "results/stats_scRNAseq/calc_cell_type_composition"


# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Load library
library(tidyverse)

# Read data
df1 <- read_csv(path_data)

# 
df1 %>%
    rename(Cell_type=types, Patient=sampleid) %>%
    group_by(Patient, Cell_type) %>%
    summarise(N_cell = n()) %>%
    group_by(Patient) %>%
    mutate(Composition = N_cell/sum(N_cell)) -> df2

# Write table
df2 %>%
    write_csv(file.path(path_outdir, "cell_type_composition.csv"))

# Draw barplot
df2 %>%
    ggplot(aes(Patient, Composition, fill=Cell_type)) +
    geom_bar(stat="identity", position = "stack") +
    theme(text=element_text(size = 16,  family="Arial")) +
    theme_bw() +
    labs(y="Cell type composition") -> g1
ggsave(file.path(path_outdir, "barplot_cell_type_composition.pdf"), g1)
ggsave(file.path(path_outdir, "barplot_cell_type_composition.svg"), g1)

# Draw barplot
df2 %>%
    ggplot(aes(Patient, N_cell, fill=Cell_type)) +
    geom_bar(stat="identity", position = "stack") +
    theme(text=element_text(size = 16,  family="Arial")) +
    theme_bw() +
    labs(y="Cell type number") -> g2
ggsave(file.path(path_outdir, "barplot_cell_type_number.pdf"), g2)
ggsave(file.path(path_outdir, "barplot_cell_type_number.svg"), g2)
