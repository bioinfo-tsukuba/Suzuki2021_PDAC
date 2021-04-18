#!/bin/bash

# Define base directory
base_dir=results/ICGC_unsupervised

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
R --vanilla -e "library(Seurat); library(tidyverse); sessionInfo()" >> ${log_file}

# Run dimensional reduction and clustering analysis
Rscript --vanilla scripts/ICGC_unsupervised/010-ICGC_unsupervised.R ${out_dir} >> ${log_file}
