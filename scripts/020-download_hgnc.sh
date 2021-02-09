#!/bin/sh

# EnsembleをGene Symbolに変換するために遺伝子データをダウンロードする

mkdir -p data/HGNC

wget -O - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt |
  sed 1d |
  awk -F "\t" 'BEGIN{OFS="\t"} {print $2,$20}' | # <- ensembleとgene symbolのみを抽出
  awk 'NF==2' |
  sort |
  join - data/HGNC/natmi_genes.txt |
cat > data/HGNC/hgnc_ensemble_symbol.txt