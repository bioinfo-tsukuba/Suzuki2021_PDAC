#!/bin/bash

# Define base directory
base_dir=results/survival_meta_analysis_plot
result=results/Fig1/LR_HR_adjPval.csv

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

# Run metap analysis
Rscript --vanilla scripts/survival_meta_analysis_plot/survival_meta_analysis_plot.R \
    ${result} ${out_dir} >> ${log_file}
