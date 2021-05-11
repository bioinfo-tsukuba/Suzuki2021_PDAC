################################################################################
# Initialization
################################################################################

library(tidyverse)
library(cowplot)

################################################################################
# Input and format
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("At least one argument is required")
}

path_result <- args[1]
out_dir <- args[2]

if(!dir.exists(out_dir)){
  dir.create(out_dir, recursive=TRUE)
}

#####################################
# Read result
df1 <- read_csv(path_result)

# Add annotation to LR pairs
df1 %>%
  mutate(`HR_PAAD-US` = case_when(
    is.na(`HR_PAAD-US`) ~ 1,
    TRUE ~ `HR_PAAD-US`
  )) %>%
  mutate(annot = case_when(
    adjPval <0.1 & `HR_PAAD-US` <1 & `HR_PACA-AU` <1 & `HR_PACA-CA` <1 ~ "0 HR<1, q<0.1",
    adjPval <0.1 & `HR_PAAD-US` >1 & `HR_PACA-AU` >1 & `HR_PACA-CA` >1 ~ "1 HR>1, q<0.1",
    TRUE ~ "2 The others"
  )) %>%
  mutate(
    HR_geometric_mean = (`HR_PAAD-US`*`HR_PACA-AU`*`HR_PACA-CA`)^(1/3),
    HR_mean = (`HR_PAAD-US`+`HR_PACA-AU`+`HR_PACA-CA`)/3
  ) -> df1

# Save csv
write_csv(df1, file.path(out_dir, "results_LR_meta.for_plot.csv"))

######################################
# Plot

df1 %>%
  group_by(annot) %>%
  summarise(n=n()) %>%
  mutate(label=sprintf("%s (%d)", annot, n)) -> df1_summary


df1 %>%
  ggplot(aes(HR_geometric_mean, -log10(adjPval), color=annot)) + 
  geom_point(size=0.5) +
  theme_half_open() +
  background_grid() +
  scale_color_manual(labels=df1_summary$label, values=c("red", "blue", "lightgray")) +
  labs(
    title="Survival meta analysis", 
    x = "Geometric mean of hazard ratios across 3 cohorts",
    y = "-log10(q-value) by Storey's method from meta p-value"
    )-> g
g

ggsave(file.path(out_dir, "geometric_mean_HR_vs_adjPval.pdf"))


df1 %>%
  ggplot(aes(HR_mean, -log10(adjPval), color=annot)) + 
  geom_point(size=0.5) +
  theme_half_open() +
  background_grid() +
  scale_color_manual(labels=df1_summary$label, values=c("red", "blue", "lightgray")) +
  labs(
    title="Survival meta analysis", 
    x = "Mean of hazard ratios across 3 cohorts",
    y = "-log10(q-value) by Storey's method from meta p-value"
    )-> g
g

ggsave(file.path(out_dir, "mean_HR_vs_adjPval.pdf"))


######################################
# sessionInfo()
sessionInfo()

