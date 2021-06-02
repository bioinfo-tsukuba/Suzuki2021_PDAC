#!/bin/sh

cat results/Fig1/LR_Pval_HR.csv |
  awk -F, 'BEGIN{OFS=","}
  $3<0.05 {
    P[$2]+=$3
    HR=$4
    if(HR<1) {
      LR[$2]++
    } else {
      LR[$2]--
    }} END {
    for(a in LR) {
      if (LR[a]==3) print a,P[a]/3,"Low"
      else if (LR[a]==-3) print a,P[a]/3,"High"
    }
  }' |
  sort -t, -k2,2n -k3,3 |
  awk -F, 'BEGIN{print "LR,mean_Pval,HR"}1' |
  cat >results/Fig1/kaplan_meier/LR_sigificance_in_all_cohorts.csv
