path_data <- "data/WeiLin_pdac10/PDAC10_meta.csv"
path_outdir <- "results/FigS3"
path_metadata <- "data/WeiLin_pdac10/Patient_data.csv"

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
  rename(Cell_type = types, Patient_ID = sampleid) %>%
  group_by(Patient_ID, Cell_type) %>%
  summarise(N_cell = n()) %>%
  group_by(Patient_ID) %>%
  mutate(Composition = N_cell / sum(N_cell), N_cell_sum = sum(N_cell)) -> df2

#
dfmeta <- read_csv(path_metadata)

#
df2 %>%
  left_join(dfmeta, by = "Patient_ID") -> df2

# Write table
df2 %>%
  write_csv(file.path(path_outdir, "cell_type_composition.csv"))

# Draw barplot
df2 %>%
  ggplot(aes(Patient_ID, Composition, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_text(size = 8, family = "Helvetica", color = "black"),
    axis.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.title = element_text(size = 6, family = "Helvetica", color = "black"),
    strip.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.key.size = unit(0.05, "inch"),
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1.5
  ) +
  labs(y = "Number of cells") -> g1
ggsave(file.path(path_outdir, "barplot_cell_type_composition.pdf"), g1,
  width = 3.5, height = 3.5
)
ggsave(file.path(path_outdir, "barplot_cell_type_composition.svg"), g1,
  width = 3.5, height = 3.5
)

# Draw barplot
df2 %>%
  ggplot(aes(Patient_ID, N_cell, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_text(size = 8, family = "Helvetica", color = "black"),
    axis.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.title = element_text(size = 6, family = "Helvetica", color = "black"),
    strip.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.key.size = unit(0.05, "inch"),
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  ) +
  labs(y = "Number of cells") -> g2
ggsave(file.path(path_outdir, "barplot_cell_type_number.pdf"), g2,
  width = 3.5, height = 3.5
)
ggsave(file.path(path_outdir, "barplot_cell_type_number.svg"), g2,
  width = 3.5, height = 3.5
)

# Draw barplot
df2 %>%
  ggplot(aes(Patient_ID, Composition, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(cols = vars(Grade), scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_text(size = 8, family = "Helvetica", color = "black"),
    axis.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.title = element_text(size = 6, family = "Helvetica", color = "black"),
    strip.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.key.size = unit(0.05, "inch"),
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 10
  ) +
  labs(y = "Cell type composition") -> g3
ggsave(file.path(path_outdir, "barplot_cell_type_composition_by_grade.pdf"), g3,
  width = 3.5, height = 3.5
)
ggsave(file.path(path_outdir, "barplot_cell_type_composition_by_grade.svg"), g3,
  width = 3.5, height = 3.5
)

# Draw barplot
df2 %>%
  ggplot(aes(Patient_ID, Composition, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(rows = vars(Cell_type), cols = vars(Grade), scales = "free", space = "free_x") +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_text(size = 8, family = "Helvetica", color = "black"),
    axis.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.title = element_text(size = 6, family = "Helvetica", color = "black"),
    strip.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.key.size = unit(0.05, "inch"),
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1.5
  ) +
  labs(y = "Cell type composition") -> g4
ggsave(file.path(path_outdir, "barplot_cell_type_composition_by_grade_and_celltype.pdf"), g4,
  width = 3.5, height = 3.5
)
ggsave(file.path(path_outdir, "barplot_cell_type_composition_by_grade_and_celltype.svg"), g4,
  width = 3.5, height = 3.5
)

# Draw barplot
df2 %>%
  ggplot(aes(Patient_ID, N_cell, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(rows = vars(Cell_type), cols = vars(Grade), scales = "free", space = "free_x") +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_text(size = 8, family = "Helvetica", color = "black"),
    axis.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.title = element_text(size = 6, family = "Helvetica", color = "black"),
    strip.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.key.size = unit(0.05, "inch"),
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1.5
  ) +
  labs(y = "Number of cells") -> g5
ggsave(file.path(path_outdir, "barplot_cell_type_number_by_grade_and_celltype.pdf"), g5,
  width = 3.5, height = 3.5
)
ggsave(file.path(path_outdir, "barplot_cell_type_number_by_grade_and_celltype.svg"), g5,
  width = 3.5, height = 3.5
)


# Draw barplot
df2 %>%
  ggplot(aes(Patient_ID, N_cell, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(cols = vars(Grade), scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_text(size = 8, family = "Helvetica", color = "black"),
    axis.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.title = element_text(size = 6, family = "Helvetica", color = "black"),
    strip.text = element_text(size = 6, family = "Helvetica", color = "black"),
    legend.key.size = unit(0.05, "inch"),
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 10
  ) +
  labs(y = "Number of cells") -> g6
ggsave(file.path(path_outdir, "barplot_cell_type_number_by_grade.pdf"), g6,
  width = 3.5, height = 3.5
)
ggsave(file.path(path_outdir, "barplot_cell_type_number_by_grade.svg"), g6,
  width = 3.5, height = 3.5
)
