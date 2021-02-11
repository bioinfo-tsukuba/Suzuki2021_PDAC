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
echo "*" > data/ICGC/.gitignore # <- ファイルが重たいのでgitignoreします

for project in PAAD-US PACA-AU PACA-CA; do
cat << FILE |
https://dcc.icgc.org/api/v1/download?fn=/current/Projects/$project/donor.$project.tsv.gz
https://dcc.icgc.org/api/v1/download?fn=/current/Projects/$project/exp_seq.$project.tsv.gz
FILE
  while read -r file; do
    wget -O - "$file" > "data/ICGC/${file##*/}"
  done
done

# 患者さんIDごとに遺伝子発現をCPM正規化します

for project in PAAD-US PACA-AU PACA-CA; do
  gzip -dc "data/ICGC/exp_seq.$project.tsv.gz" |
    grep -v "icgc_donor_id" |
    awk -F "\t" 'BEGIN {OFS="\t"} {
      id[NR]=$1; gene[NR]=$8; value[NR]=$10; sum[$1]+=$10
      } END {
      for(nr in id) {
        key=id[nr]
        cpm=value[nr]/sum[key]*1000000
        print gene[nr], key, cpm
      }
      }' > data/ICGC/exp_seq_"$project".txt
done
