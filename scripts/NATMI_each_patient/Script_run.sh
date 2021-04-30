#!/bin/bash

commit=`git log --format="%H" -n 1 | cut -c1-7`
date=$(date '+%Y%m%d-%H%M%S')

# Run R for making input data for each patient
mkdir -p processed_data/NATMI_each_patient
echo "git commit: ${commit}" > processed_data/NATMI_each_patient/log.txt
echo "git diff --cached:" >> processed_data/NATMI_each_patient/log.txt
git diff --cached >> processed_data/NATMI_each_patient/log.txt

Rscript --vanilla scripts/NATMI_each_patient/010-make_NATMI_input.R >> processed_data/NATMI_each_patient/log.txt

# Define base directory
base_dir=results/NATMI_each_patient

# Set output directory
out_dir=${base_dir}
log_file=${out_dir}/log.txt

# Make output directory
mkdir -p $out_dir

# Log
echo ${out_dir} > ${log_file}
echo "git commit: ${commit}" >> ${log_file}
echo "git diff --cached:" >> ${log_file}
git diff --cached >> ${log_file}

# Run 020-run_NATMI_each_patient.sh
bash scripts/NATMI_each_patient/020-run_NATMI_each_patient.sh 2>&1 >> results/NATMI_each_patient/log.txt
