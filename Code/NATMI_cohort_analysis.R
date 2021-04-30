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

LR_common_pval <- NATMI_LR_common %>% dplyr::select(c("lr_pair", "pval_US", "pval_CA", "pval_AU" ))

write.csv(LR_common_pval, "/Users/sayakasuzuki/Desktop/SSD/results/NATMI_LR_common_pval.csv")

#read NATMI_meta
NATMI_meta <- read.csv("/Users/sayakasuzuki/Desktop/SSD/results/survival_meta_analysis/20210417-115137-ab27039/results_LR_meta.csv")

#filter hr<1,qpval<0.1
NATMI_meta %>% filter(hr_PAAD.US<1, hr_PACA.AU<1, hr_PACA.CA<1, meta_pval<0.1) ->NATMI_meta_filtered_pval

write.csv(NATMI_meta_filtered_pval, "/Users/sayakasuzuki/Desktop/SSD/results/survival_meta_analysis//NATMI_meta_filtered_pval.csv")

#filter BH<0.1, hr<1
NATMI_meta %>% filter(meta_qval_BH < 0.1, hr_PAAD.US < 1, hr_PACA.AU< 1, hr_PACA.CA< 1) ->NATMI_meta_filtered_BH
write.csv(NATMI_meta_filtered_BH, "/Users/sayakasuzuki/Desktop/SSD/results/survival_meta_analysis//NATMI_meta_filtered_BH.csv")

#filter Holm<0.1, hr<1
NATMI_meta %>% filter(meta_qval_holm < 0.1, hr_PAAD.US < 1, hr_PACA.AU< 1, hr_PACA.CA< 1) ->NATMI_meta_filtered_Holm
write.csv(NATMI_meta_filtered_Holm, "/Users/sayakasuzuki/Desktop/SSD/results/survival_meta_analysis//NATMI_meta_filtered_Holm.csv")

#filter Storey<0.1, hr<1
NATMI_meta %>% filter(meta_qval_storey< 0.1, hr_PAAD.US < 1, hr_PACA.AU< 1, hr_PACA.CA< 1) ->NATMI_meta_filtered_Storey
write.csv(NATMI_meta_filtered_Storey, "/Users/sayakasuzuki/Desktop/SSD/results/survival_meta_analysis//NATMI_meta_filtered_Storey.csv")

#join 3 filtered dataframe
df1 <-dplyr::inner_join(NATMI_meta_filtered_Storey,NATMI_meta_filtered_BH, by="LR")
NATMI_meta_common <- dplyr::inner_join(df1, NATMI_meta_filtered_Holm, by="LR")

#joinで重複した列を削除
NATMI_meta_common <-NATMI_meta_common[!duplicated(as.list(NATMI_meta_common))]

write.csv(NATMI_meta_common, "/Users/sayakasuzuki/Desktop/SSD/results/survival_meta_analysis//NATMI_meta_common.csv")





