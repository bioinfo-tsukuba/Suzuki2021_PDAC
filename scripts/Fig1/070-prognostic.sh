#!/bin/sh

mkdir -p results/Fig1/prognostic/

sed 1d data/prognostic_genes/favorable.tsv |
  cut -f 1 |
  sort -u >tmp_favorable.txt

sed 1d data/prognostic_genes/unfavorable.tsv |
  cut -f 1 |
  sort -u >tmp_unfavorable.txt

sed 1d results/Fig1/LR_adjPval_meanHR_screened.csv |
  awk -F, '$NF<1 {sub("->","\n",$1); print $1}' |
  sort -u |
  join - tmp_favorable.txt |
  cat >results/Fig1/prognostic/favorable_genes.txt

sed 1d results/Fig1/LR_adjPval_meanHR_screened.csv |
  awk -F, '$NF>1 {sub("->","\n",$1); print $1}' |
  sort -u |
  join - tmp_unfavorable.txt |
  cat >results/Fig1/prognostic/unfavorable_genes.txt

rm tmp_

# # joined favorable genes were really small... 予後が良いLR遺伝子が少ないからかな？

# cat data/NATMI_LR.csv |
#   awk '{sub("->","\n",$1); print $1}' |
#   sort -u |
#   join - tmp_favorable.txt |
#   wc -l

# cat data/NATMI_LR.csv |
#   awk '{sub("->","\n",$1); print $1}' |
#   sort -u |
#   join - tmp_unfavorable.txt |
#   wc -l
