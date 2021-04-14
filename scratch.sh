#!/bin/sh

cat results_LR.csv |
  awk -F, '$3<1' >tmp.csv

for project in PAAD-US PACA-AU PACA-CA; do
  Rscript library/CCI_pval.R \
    data/ICGC/survival_"$project".csv.gz \
    tests/MD3/data/CCI.csv |
    sed "s/^/${project},/"
done >results_CCI.csv

cat results_CCI.csv | grep "CAF->ETC"
awk -F, '$3<0.05 && $4 < 1' |
  cut -d, -f2 |
  sort | uniq -c
