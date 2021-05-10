#!/bin/sh

mkdir -p results/Fig1

# Donwload and format NATMI LR pairs
wget -O - https://raw.githubusercontent.com/asrhou/NATMI/master/lrdbs/lrc2p.csv |
  tr "\r" "\n" |
  sed 1d |
  sed "s/,/->/" |
  grep ^ >data/NATMI_LR.csv

# Calculate Pvalue and HR
find data/ICGC/survival_*.csv.gz |
  while read -r line; do
    cohort="$(echo ${line##*_} | sed "s/.csv.gz//")"
    Rscript ./library/LR_Pval_HR.R \
      "$line" \
      data/NATMI_LR.csv |
      sed 1d |
      sed "s/^/${cohort},/"
  done |
  awk 'BEGIN{print "Cohort,LR,Pval,HR"}1' >results/Fig1/LR_Pval_HR.csv
