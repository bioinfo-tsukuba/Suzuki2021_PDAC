#!/bin/sh

# データソースは以下の3プロジェクトです
#（https://docs.icgc.org/submission/projects/）
# ====================================
# PAAD-US	Pancreatic Cancer - TCGA, US	US
# PACA-AU	Pancreatic Cancer Endocrine Neoplasms- AU	Australia
# PACA-CA	Pancreatic Cancer - CA	Canada
# ====================================
# *PACA-CNは遺伝子発現データがないので解析対象から除外しました

mkdir -p data/ICGC

for project in PAAD-US PACA-AU PACA-CA; do

cat << FILE |
https://dcc.icgc.org/api/v1/download?fn=/current/Projects/$project/donor.$project.tsv.gz
https://dcc.icgc.org/api/v1/download?fn=/current/Projects/$project/exp_seq.$project.tsv.gz
FILE
  while read -r file; do
    wget -O - "$file" > "${file##*/}"
    mv "${file##*/}" "data/ICGC"
  done

done
