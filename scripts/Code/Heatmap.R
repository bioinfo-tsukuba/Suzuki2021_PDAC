# Downloding Complex Heatmap
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("ComplexHeatmap")


library(tidyverse)
library(ggplot2)
library(Cairo)


# read csv
df_summary <- read_csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_each_patient/summary_number_of_patient.csv")

ggplot(df_summary, aes(x = LR, y = cell_type_pair, fill = number_of_patient)) +
  geom_tile() -> g
g <- g + theme_bw()

g <- g + theme(
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  strip.background = element_rect(fill = "white", colour = "white"),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
)

g <- g + scale_fill_gradientn("Number of Patient", colours = rev(brewer.pal(9, "Spectral")), na.value = "white")
