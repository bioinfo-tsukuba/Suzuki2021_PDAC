path_data1 <- "results/stats_scRNAseq/calc_cell_type_composition/cell_type_composition.csv"
path_data2 <- "results/NATMI_each_patient/Dataframe_AtoE/Dataframe_D.csv"
path_outdir <- "results/NATMI_each_patient/compare_with_cell_type_composition"

# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Load library
library(tidyverse)

# Read data
df1 <- read_csv(path_data1)
df2 <- read_csv(path_data2)

#
df2 %>%
  separate(cell_type_pair, sep = "->", into = c("From", "To"), remove = FALSE) %>%
  left_join(df1, by = c("From" = "Cell_type", "Patient")) %>%
  left_join(df1, by = c("To" = "Cell_type", "Patient")) -> df3
df3 %>%
  rename(N_cell_sending = N_cell.x, N_cell_receiving = N_cell.y) -> df3

# Write table
df3 %>%
  write_csv(file.path(path_outdir, "NATMI_and_cell_type_number.csv"))


# df3
df3 %>%
  ggplot(aes(log10(N_cell_sending * N_cell_receiving), NATMI_LR_ALL, col = Patient)) +
  geom_point() +
  theme(text = element_text(size = 16)) +
  theme_bw() +
  labs(
    x = "log10(N_sending * N_receiving)",
    y = "Number of LR pairs detected by NATMI",
    title = "Detected LR pairs for each cell type pair and each patient"
  ) -> g1
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_ALL.pdf"), g1)
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_ALL.svg"), g1)

df3 %>%
  ggplot(aes(log10(N_cell_sending * N_cell_receiving), HRL, col = Patient)) +
  geom_point() +
  theme(text = element_text(size = 16)) +
  theme_bw() +
  labs(
    x = "log10(N_sending * N_receiving)",
    y = "Number of LR pairs detected by NATMI",
    title = "Detected LR pairs for each cell type pair and each patient"
  ) -> g2
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_HRL.pdf"), g2)
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_HRL.svg"), g2)

df3 %>%
  ggplot(aes(log10(N_cell_sending * N_cell_receiving), HRH, col = Patient)) +
  geom_point() +
  theme(text = element_text(size = 16)) +
  theme_bw() +
  labs(
    x = "log10(N_sending * N_receiving)",
    y = "Number of LR pairs detected by NATMI",
    title = "Detected LR pairs for each cell type pair and each patient"
  ) -> g3
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_HRH.pdf"), g3)
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_HRH.svg"), g3)
