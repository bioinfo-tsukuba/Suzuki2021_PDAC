################################################################################
# Initialization
################################################################################

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, metap, qvalue)

################################################################################
# Input and format
################################################################################

df <- read_csv("results/Fig1/LR_Pval_HR.csv")

# Convert to wider data
df_wider <- df %>%
  pivot_wider(id_cols = LR, names_from = Cohort, values_from = c(Pval, HR)) %>%
  na.omit()

# Calculate adjusted p-values
mat_pval <- df_wider %>%
  select(starts_with("Pval")) %>%
  as.matrix()

meta_pval <- apply(mat_pval, 1, function(x) sump(x)$p)
meta_pval_storey <- qvalue(meta_pval)$qvalues

df_meta_pval <- df_wider %>%
  mutate(meta_Pval = meta_pval) %>%
  mutate(adjPval = meta_pval_storey) %>%
  inner_join(df) %>%
  select(LR, HR, Pval, meta_Pval, adjPval)

# Filter LR by p-val and HR
df_result <- df_meta_pval %>%
  # Pval after metaP and Storey < 0.1
  filter(adjPval < 0.1) %>%
  group_by(LR) %>%
  mutate(high_HR = sum(HR > 1)) %>%
  mutate(low_HR = sum(HR < 1)) %>%
  # Extract HR>1 or HR<1 in every three cohorts
  filter(high_HR == 3 || low_HR == 3) %>%
  group_by(LR) %>%
  mutate(meanHR = mean(HR)) %>%
  select(LR, adjPval, meanHR) %>%
  distinct() %>%
  arrange(adjPval)

################################################################################
# Output
################################################################################

write_csv(df_result, "results/Fig1/LR_adjPval_meanHR_screened.csv")
