################################################################################
# Initialization
################################################################################

library(metap)
library(tidyverse)
library(qvalue)

################################################################################
# Input and format
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("At least one argument is required")
}

#####################################
# Read ICGC data
lapply(args, function(path_survival_result) {
  read_csv(path_survival_result) %>%
    mutate(cohort = gsub(".csv", "", gsub("results_LR_", "", basename(path_survival_result))))
}) %>%
  bind_rows() -> df1

# Convert to wider data
df1 %>%
  pivot_wider(
    id_cols = LR,
    names_from = cohort,
    values_from = c(pval, hr, median_day_low, median_day_high, diff_day)
  ) -> df2

# Calculate p-values
df2 %>%
  select(starts_with("pval")) %>%
  as.matrix() -> mat_pval

meta_pval <- numeric(nrow(mat_pval))

for (i in 1:nrow(mat_pval)) {
  if (all(!is.na(mat_pval[i, ]))) {
    meta_pval[i] <- sump(mat_pval[i, ])$p
  } else {
    # If there is an NA value, the meta p-value is set to 1.
    meta_pval[i] <- 1.0
  }
}

df2 %>%
  mutate(meta_pval = meta_pval) -> df3

# Calculate FDR for correction for multiple testing
df3 %>%
  mutate(
    meta_qval_BH = p.adjust(meta_pval, method = "BH"),
    meta_qval_holm = p.adjust(meta_pval, method = "holm"),
    meta_qval_storey = qvalue(df3$meta_pval)$qvalues
  ) -> df4


################################################################################
# Output
################################################################################

cat(format_csv(df4))
