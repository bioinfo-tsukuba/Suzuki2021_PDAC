#!/bin/sh

# 患者さんID, 性別, 生死, 観測期間, 遺伝子名, 遺伝子発現の6コラムのCSVファイルをつくります

mkdir -p data/ICGC

gzip -dc data/ICGC/donor.* |
  grep -v "icgc_donor_id" |
  cut -f 1,5,6,17,18 |
  awk 'NF==5' |
  sort > tmp_donor

gzip -dc data/ICGC/exp_seq.* |
  grep -v "icgc_donor_id" |
  awk -F "\t" 'BEGIN{OFS="\t"}{print $8,$1,$9}' |
  sort |
  join -a 1 - data/HGNC/natmi_symbol.txt |
  join -a 1 - data/HGNC/natmi_ensemble.txt |
  awk '/LR$/ {print $1,$2,$3}' |
  join -a 1 - data/HGNC/natmi_ensemble_symbol.txt |
  awk 'NF==4 {$1=$4} {print $2,$1,$3}' |
  sort -t " " |
  join tmp_donor - |
  awk 'BEGIN {OFS=","; print "id", "sex", "status", "time", "gene", "exp"}
    {print $1,$2,$3,$4,$5,$6}' |
  gzip -c > data/ICGC/survival.csv.gz

rm tmp_donor
