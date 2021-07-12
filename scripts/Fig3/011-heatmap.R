
################################################################################
# Initialization
################################################################################
options(repos = "https//cran.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pacman::p_load(tidyverse, ComplexHeatmap, circlize)
system("mkdir -p results/Fig3/")


################################################################################
# Import and format data
################################################################################
df <- read_csv("results/NATMI_each_patient/summary_number_of_patient.csv")

# Convert tidy data into matrix
mat1 <- df %>%
  pivot_wider(
    names_from = cell_type_pair,
    values_from = number_of_patient,
    values_fill = 0
  ) %>%
  select(-LR) %>%
  as.matrix()

rownames(mat1) <- df$LR %>% unique()

################################################################################
# Main
################################################################################

pdf("results/Fig3/heatmap_by_patient_number.pdf")
f1 <- colorRamp2(seq(min(mat1), max(mat1), length = 2), c("#EEEEEE", "blue"))
# heatmap
Heatmap(mat1,
  col = f1,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5)
)
dev.off()