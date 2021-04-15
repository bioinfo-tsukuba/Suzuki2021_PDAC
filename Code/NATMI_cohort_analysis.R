library(tidyverse)
library(ggplot2)

US <- read.csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_lrpair_result_US.csv")
CA <- read.csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_lrpair_result_CA.csv")
AU <- read.csv("/Users/sayakasuzuki/Desktop/SSD/results/NATMI_lrpair_result_AU.csv")

top20_US <- US %>%
  top_n(-20, pval_US)

top20_CA <- CA %>%
  top_n(-20, pval_CA)

top20_AU <- AU %>%
  top_n(-20, pval_AU)

#join 3 cohorts
LR_common <- dplyr::inner_join(US,CA, by="lr_pair")

NATMI_LR_common <- dplyr::inner_join(LR_common, AU, by="lr_pair")

write.csv(NATMI_LR_common, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_LR_common.csv")