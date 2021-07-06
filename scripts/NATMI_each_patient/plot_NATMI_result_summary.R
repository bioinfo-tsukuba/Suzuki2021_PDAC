
# Define path
in_file <- "results/NATMI_each_patient/cell_type_sumamry.csv"
path_outdir <- "results/NATMI_each_patient/plot_NATMI_result_summary"

# Load libraries
library(tidyverse)

# Make output directory
if (!dir.exists(path_outdir)) {
  dir.create(path_outdir, recursive = TRUE)
}

# Read data
df1 <- read_csv(in_file)
df1_long <- df1 %>%
    pivot_longer(P01_amount:P10_amount, names_to="Sample_ID", values_to="LR_pair") %>%
    mutate(Sample_ID=gsub("_amount", "", Sample_ID))

# Get highest number
df1_long %>%
    group_by(cell_type_pair) %>%
    summarise(sum_LR_pair = sum(LR_pair, na.rm=TRUE)) %>%
    ungroup() %>%
    summarise(max_sum_LR_pair=max(sum_LR_pair, na.rm=TRUE)) %>%
    { .$max_sum_LR_pair} -> max_sum_LR_pair

# Plot bar plot
df1_long %>%
    ggplot(aes(x=cell_type_pair, y=LR_pair, fill=Sample_ID)) +
    geom_bar(stat="identity", position = "stack") +
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
        aspect.ratio = 0.3
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, max_sum_LR_pair*1.05)) + # max_sum_LR_pair
    theme(legend.position="bottom") +
    labs(x="Cell type pair", y="") -> g1

ggsave(file.path(path_outdir, "barplot_NATMI_LRpairs.svg"), g1,
    width = 4, height = 4
)
ggsave(file.path(path_outdir, "barplot_NATMI_LRpairs.pdf"), g1,
    width = 4, height = 4
)

# Session info
sessionInfo()
