#!/bin/bash

# Define base directory
base_dir=results/survival_meta_analysis

# Set output directory
commit=`git log --format="%H" -n 1 | cut -c1-7`
date=$(date '+%Y%m%d-%H%M%S')
out_dir=${base_dir}/${date}-${commit}
log_file=${out_dir}/log.txt

# Make output directory
mkdir -p $out_dir

# Log
echo ${out_dir} > ${log_file}
echo "git commit: ${commit}" >> ${log_file}
echo "git diff --cached:" >> ${log_file}
git diff --cached >> ${log_file}
R --vanilla -e "pacman::p_load(survival, survminer, broom, tidyverse); sessionInfo()" >> ${log_file}
R --vanilla -e "library(metap); library(tidyverse); library(qvalue); sessionInfo()" >> ${log_file}

# Run survival analysis
for cohort in "PAAD-US" "PACA-AU" "PACA-CA"
do
    Rscript --vanilla ./library/LR_pval.R \
    data/ICGC/survival_${cohort}.csv.gz \
    tests/MD3/data/LR.csv > ${out_dir}/results_LR_${cohort}.csv
done

# Run metap analysis
results=`ls ${out_dir}/results_LR_*.csv`
Rscript --vanilla scripts/survival_meta_analysis/metap.R \
    ${results} > ${out_dir}/results_LR_meta.csv
