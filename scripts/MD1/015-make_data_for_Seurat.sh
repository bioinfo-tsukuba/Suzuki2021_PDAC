#!/bin/bash

# Make output directory
mkdir -p data/WeiLin_pdac10/For_Seurat

# Make features.tsv.gz
gzip -c data/WeiLin_pdac10/genes.tsv >data/WeiLin_pdac10/For_Seurat/features.tsv.gz
cat data/WeiLin_pdac10/genes.tsv | awk 'BEGIN{OFS="\t"}{print $1,$1}' | gzip -c >data/WeiLin_pdac10/For_Seurat/features.tsv.gz

# Make barcodes.tsv.gz
awk -F"," 'NR>1{gsub(/\"/, ""); print $1}' data/WeiLin_pdac10/PDAC10_meta.csv | gzip -c >data/WeiLin_pdac10/For_Seurat/barcodes.tsv.gz

# Make matrix.mtx.gz
cp data/WeiLin_pdac10/PDAC10.mtx.gz data/WeiLin_pdac10/For_Seurat/matrix.mtx.gz
