#
path_data1 <- "results/stats_scRNAseq/calc_cell_type_composition/cell_type_composition.csv"
path_data2 <- "results/NATMI_each_patient/Dataframe_AtoE/Dataframe_D.csv"
#
path_LR <- "results/Fig1/LR_adjPval_meanHR_screened.csv"
path_LR_all <- "results/Fig1/LR_HR_adjPval.csv"
#
path_outdir <- "results/NATMI_each_patient/compare_with_cell_type_composition"
#
path_metadata <- "data/WeiLin_pdac10/Patient_data.csv"


# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Load library
library(tidyverse)
library(ggrepel)


# Read data
df1 <- read_csv(path_data1)
df2 <- read_csv(path_data2)

# Read metadata
dfmeta <- read_csv(path_metadata) %>%
  filter(Primary_or_Metasitasis == "Primary")

#
#
df2 %>%
  separate(cell_type_pair, sep = "->", into = c("From", "To"), remove = FALSE) %>%
  left_join(df1, by = c("From" = "Cell_type", "Patient" = "Patient_ID")) %>%
  left_join(df1, by = c("To" = "Cell_type", "Patient" = "Patient_ID")) -> df3
df3 %>%
  rename(N_cell_sending = N_cell.x, N_cell_receiving = N_cell.y) -> df3

# Write table
df3 %>%
  write_csv(file.path(path_outdir, "NATMI_and_cell_type_number.csv"))


# df3
df3 %>%
  ggplot(aes(log10(N_cell_sending * N_cell_receiving), NATMI_LR_ALL, col = Patient)) +
  geom_point(size = 0.5) +
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
        aspect.ratio = 0.7
  ) +
  labs(
    x = "log10(N_sending * N_receiving)",
    y = "Number of LR pairs detected by NATMI",
    caption = "Detected LR pairs for each cell type pair and each patient"
  ) -> g1
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_ALL.pdf"), g1, 
    width = 4, height = 4)
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_ALL.svg"), g1,
    width = 4, height = 4)

df3 %>%
  ggplot(aes(log10(N_cell_sending * N_cell_receiving), HRL, col = Patient)) +
  geom_point(size = 0.5) +
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
        aspect.ratio = 0.7
  ) +
  labs(
    x = "log10(N_sending * N_receiving)",
    y = "Number of LR pairs detected by NATMI",
    caption = "Detected LR pairs for each cell type pair and each patient"
  ) -> g2
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_HRL.pdf"), g2,
    width = 4, height = 4)
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_HRL.svg"), g2,
    width = 4, height = 4)

df3 %>%
  ggplot(aes(log10(N_cell_sending * N_cell_receiving), HRH, col = Patient)) +
  geom_point(size = 0.5) +
  geom_point(size = 0.5) +
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
        aspect.ratio = 0.7
  ) +
  labs(
    x = "log10(N_sending * N_receiving)",
    y = "Number of LR pairs detected by NATMI",
    caption = "Detected LR pairs for each cell type pair and each patient"
  ) -> g3
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_HRH.pdf"), g3,
    width = 4, height = 4)
ggsave(file.path(path_outdir, "scatter_NATMI_and_cell_type_number_HRH.svg"), g3,
    width = 4, height = 4)

######################################
# Count the number of LR pairs
dflr <- read_csv(path_LR)
dflr %>%
  filter(meanHR < 1) %>%
  nrow() -> n_HRL
dflr %>%
  filter(meanHR > 1) %>%
  nrow() -> n_HRH

# Count the number of all LR pairs
dflrall <- read_csv(path_LR_all)
dflrall %>%
  nrow() -> n_LR_all

#
df3 %>%
  group_by(cell_type_pair) %>%
  summarise(
    mean_NormHRL = mean(NormHRL),
    mean_NormHRH = mean(NormHRH),
    adjusted_mean_NormHRL = mean(NormHRL) / n_HRL * n_LR_all,
    adjusted_mean_NormHRH = mean(NormHRH) / n_HRH * n_LR_all,
    adjusted_mean_NormHRL2 = mean(NormHRL / log10((N_cell_sending + 1) * (N_cell_receiving + 1)) * log10(2)) / n_HRL * n_LR_all,
    adjusted_mean_NormHRH2 = mean(NormHRH / log10((N_cell_sending + 1) * (N_cell_receiving + 1)) * log10(2)) / n_HRH * n_LR_all,
    adjusted_mean_NormHRL3 = mean(HRL / log10((N_cell_sending + 1) * (N_cell_receiving + 1)) * log10(2)) / n_HRL,
    adjusted_mean_NormHRH3 = mean(HRH / log10((N_cell_sending + 1) * (N_cell_receiving + 1)) * log10(2)) / n_HRH,
    cell_number_product_log10 = mean(log10((N_cell_sending + 1) * (N_cell_receiving + 1)))
  ) -> df4



df4 %>%
  ggplot(aes(adjusted_mean_NormHRL, adjusted_mean_NormHRL2)) +
  geom_point(alpha = 0.3)

df4 %>%
  ggplot(aes(adjusted_mean_NormHRH, adjusted_mean_NormHRH2)) +
  geom_point(alpha = 0.3)

df3 %>%
  ggplot(aes(log10(N_cell_sending * N_cell_receiving), NormHRL / n_HRL * n_LR_all)) +
  geom_point()
df3 %>%
  ggplot(aes(log10(N_cell_sending * N_cell_receiving), NormHRH / n_HRH * n_LR_all)) +
  geom_point()

################
# Scatter plot
min_value <- min(c(df4$adjusted_mean_NormHRL2, df4$adjusted_mean_NormHRH2))
max_value <- max(c(df4$adjusted_mean_NormHRL2, df4$adjusted_mean_NormHRH2))
df4 %>%
  ggplot(aes(adjusted_mean_NormHRL2, adjusted_mean_NormHRH2, label = cell_type_pair)) +
  geom_point(aes(color = adjusted_mean_NormHRL2 / adjusted_mean_NormHRH2 > 2)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0, slope = 2, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1 / 2, linetype = "dashed") +
  geom_text_repel(
    data = df4 %>% filter(adjusted_mean_NormHRH2 / adjusted_mean_NormHRL2 > 2 | adjusted_mean_NormHRL2 / adjusted_mean_NormHRH2 > 2)
  ) +
  theme(legend.position = "bottom") +
  lims(x = c(min_value, max_value), y = c(min_value, max_value)) +
  labs(
    title = "Adjusted mean enrichment of LR pairs",
    x = "LR pairs with HR<1", y = "LR pairs with HR>1"
  ) -> g1
g1

################
# Scatter plot
min_value <- min(c(df4$adjusted_mean_NormHRL3, df4$adjusted_mean_NormHRH3))
max_value <- max(c(df4$adjusted_mean_NormHRL3, df4$adjusted_mean_NormHRH3))
df4 %>%
  ggplot(aes(adjusted_mean_NormHRL3, adjusted_mean_NormHRH3, label = cell_type_pair)) +
  geom_point(aes(color = adjusted_mean_NormHRL2 / adjusted_mean_NormHRH2 > 2)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0, slope = 2, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1 / 2, linetype = "dashed") +
  geom_text_repel(
    data = df4 %>% filter(adjusted_mean_NormHRH3 / adjusted_mean_NormHRL3 > 2 | adjusted_mean_NormHRL3 / adjusted_mean_NormHRH3 > 2)
  ) +
  theme(legend.position = "bottom") +
  lims(x = c(min_value, max_value), y = c(min_value, max_value)) +
  labs(
    title = "Adjusted mean enrichment of LR pairs",
    x = "LR pairs with HR<1", y = "LR pairs with HR>1"
  ) -> g1
g1


###################
df4 %>%
  ggplot(aes(reorder(cell_type_pair, adjusted_mean_NormHRL3), adjusted_mean_NormHRL3)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(text = element_text(size = 16)) +
  theme_bw() +
  labs(x = "Cell type pair")

df4 %>%
  ggplot(aes(reorder(cell_type_pair, adjusted_mean_NormHRL2), adjusted_mean_NormHRL2)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(text = element_text(size = 16)) +
  theme_bw() +
  labs(x = "Cell type pair")

df4 %>%
  ggplot(aes(reorder(cell_type_pair, adjusted_mean_NormHRL), adjusted_mean_NormHRL)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(text = element_text(size = 16)) +
  theme_bw() +
  labs(x = "Cell type pair")

#########################
df4 %>%
  separate(cell_type_pair, sep = "->", into = c("Sending", "Receiving"), remove = FALSE) %>%
  ggplot(aes(Sending, Receiving, fill = adjusted_mean_NormHRL)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "blue")

df4 %>%
  separate(cell_type_pair, sep = "->", into = c("Sending", "Receiving"), remove = FALSE) %>%
  ggplot(aes(Sending, Receiving, fill = adjusted_mean_NormHRL2)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "blue")

df4 %>%
  separate(cell_type_pair, sep = "->", into = c("Sending", "Receiving"), remove = FALSE) %>%
  ggplot(aes(Sending, Receiving, fill = adjusted_mean_NormHRL3)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "blue")

df4 %>%
  separate(cell_type_pair, sep = "->", into = c("Sending", "Receiving"), remove = FALSE) %>%
  ggplot(aes(Sending, Receiving, fill = cell_number_product_log10)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "blue")

#############################################################################################
#############################################################################################
#############################################################################################

#
df3 %>%
  left_join(dfmeta, by = c("Patient" = "Patient_ID")) %>%
  group_by(cell_type_pair, Grade) %>%
  summarise(
    mean_NormHRL = mean(NormHRL),
    mean_NormHRH = mean(NormHRH),
    adjusted_mean_NormHRL = mean(NormHRL) / n_HRL * n_LR_all,
    adjusted_mean_NormHRH = mean(NormHRH) / n_HRH * n_LR_all,
    adjusted_mean_NormHRL2 = mean(NormHRL / log10((N_cell_sending + 1) * (N_cell_receiving + 1)) * log10(2)) / n_HRL * n_LR_all,
    adjusted_mean_NormHRH2 = mean(NormHRH / log10((N_cell_sending + 1) * (N_cell_receiving + 1)) * log10(2)) / n_HRH * n_LR_all,
    adjusted_mean_NormHRL3 = mean(HRL / log10((N_cell_sending + 1) * (N_cell_receiving + 1)) * log10(2)) / n_HRL,
    adjusted_mean_NormHRH3 = mean(HRH / log10((N_cell_sending + 1) * (N_cell_receiving + 1)) * log10(2)) / n_HRH,
    cell_number_product_log10 = mean(log10((N_cell_sending + 1) * (N_cell_receiving + 1)))
  ) -> df5

df5 %>%
  separate(cell_type_pair, sep = "->", into = c("Sending", "Receiving"), remove = FALSE) %>%
  ggplot(aes(Sending, Receiving, fill = cell_number_product_log10)) +
  geom_tile() +
  facet_wrap(~Grade) +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "blue")
df5 %>%
  separate(cell_type_pair, sep = "->", into = c("Sending", "Receiving"), remove = FALSE) %>%
  ggplot(aes(Sending, Receiving, fill = adjusted_mean_NormHRL)) +
  geom_tile() +
  facet_wrap(~Grade) +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "blue")
df5 %>%
  separate(cell_type_pair, sep = "->", into = c("Sending", "Receiving"), remove = FALSE) %>%
  ggplot(aes(Sending, Receiving, fill = adjusted_mean_NormHRL2)) +
  geom_tile() +
  facet_wrap(~Grade) +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "blue")
df5 %>%
  separate(cell_type_pair, sep = "->", into = c("Sending", "Receiving"), remove = FALSE) %>%
  ggplot(aes(Sending, Receiving, fill = adjusted_mean_NormHRL3)) +
  geom_tile() +
  facet_wrap(~Grade) +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "blue")

df5 %>%
  separate(cell_type_pair, sep = "->", into = c("Sending", "Receiving"), remove = FALSE) %>%
  ggplot(aes(Sending, Receiving, fill = adjusted_mean_NormHRL / cell_number_product_log10)) +
  geom_tile() +
  facet_wrap(~Grade) +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "blue")
