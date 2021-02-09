#!/bin/sh

# 患者さんID, 性別, 生死, 観測期間, 遺伝子名, 遺伝子発現の6コラムのCSVファイルをつくります

mkdir -p data/ICGC
echo "*" > data/ICGC/.gitignore # <- ファイルが重たいのでgitignoreします

gzip -dc data/ICGC/donor.* |
  grep -v "icgc_donor_id" |
  cut -f 1,5,6,17,18 |
  awk 'NF==5' |
  sort |
cat > tmp_donor

gzip -dc data/ICGC/exp_seq.* |
  grep -v "icgc_donor_id" |
  cut -f 1,8,9 |
  sort |
  join tmp_donor - > tmp_
  awk 'BEGIN{OFS=","; print "id", "sex", "status", "time", "gene", "exp"}
  {print $1,$2,$3,$4,$5,$6}' |
gzip -c > data/ICGC/survival.csv.gz

rm tmp_donor
