#!/bin/sh

mkdir -p data/HGNC

wget -O - https://asrhou.github.io/NATMI/ > tmp_natmi_genes.html

cat tmp_natmi_genes.html |
  awk 'BEGIN{RS="<td class=\"col1\">"} {$1=$1}1' |
  awk '{sub("<td class=\"col4\"> ", "\n")}1' |
  cut -d " " -f 1 |
  sort -u > data/HGNC/natmi_genes.txt

rm tmp_natmi_genes.html