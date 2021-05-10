#!/bin/sh

mkdir -p results/Fig1_ICGC/kaplan_meier/

cat results/Fig1_ICGC/LR_selected.csv |
  sed 1d |
  awk -F, '$NF=="HR>1" {print $1}' >tmp_high

cat results/Fig1_ICGC/LR_selected.csv |
  sed 1d |
  awk -F, '$NF=="HR<1" {print $1}' >tmp_low

find data/ICGC/survival_*.csv.gz |
  while read -r line; do
    cohort="$(echo ${line##*_} | sed "s/.csv.gz//")"
    Rscript ./library/LR_plot.R "$line" tmp_high "results/Fig1_ICGC/kaplan_meier/${cohort}_highHR.pdf"
    Rscript ./library/LR_plot.R "$line" tmp_low "results/Fig1_ICGC//kaplan_meier/${cohort}_lowHR.pdf"
  done

rm tmp*
