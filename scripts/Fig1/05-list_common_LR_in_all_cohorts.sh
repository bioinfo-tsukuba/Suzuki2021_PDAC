#!/bin/sh

cat results/Fig1/LR_Pval_HR.csv |
  awk -F, 'BEGIN{OFS=","}
  $3<0.05 {
    HR=$NF
    if(HR<1) {
      LR[$2]++
    } else {
      LR[$2]--
    }} END {
    for(a in LR) {
      if (LR[a]==3) print a,"HR<1"
      else if (LR[a]==-3) print a,"HR>1"
    }
  }' |
  sort -t, -k2,2 |
  cat >results/Fig1/kaplan_meier/LR_sigificance_in_all_cohorts.csv
