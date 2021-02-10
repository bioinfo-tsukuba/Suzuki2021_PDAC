#!/bin/sh

# EnsembleをGene Symbolに変換するために遺伝子データをダウンロードする

mkdir -p data/HGNC

wget -O - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt |
  sed 1d |
  cut -f 2,20 | # <- ensembleとgene symbolのみを抽出
  awk 'NF==2' |
  sort |
  join - data/HGNC/natmi_genes.txt |
  awk '{print $2,$1}' |
  sort -t " " > data/HGNC/natmi_ensemble_symbol.txt

cut -d " " -f 1 data/HGNC/natmi_ensemble_symbol.txt | sort | awk '{print $1,"LR"}' > data/HGNC/natmi_ensemble.txt
cut -d " " -f 2 data/HGNC/natmi_ensemble_symbol.txt | sort | awk '{print $1,"LR"}' > data/HGNC/natmi_symbol.txt
