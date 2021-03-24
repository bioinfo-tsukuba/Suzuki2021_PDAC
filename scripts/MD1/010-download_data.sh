#!/bin/sh

# Dependencies; wget, gzip

## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154778
for GSM in $(awk 'BEGIN{for(i=4679532;i<=4679541;i++) print i}'); do
  echo GSM"$GSM"...
  wget -qO - "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM$GSM" |
    grep ftp |
    grep "$GSM" |
    sed "s/^.*href=\"ftp/\"ftp/" |
    sed "s/\">.*/\"/" |
    sed "s/^/wget -q /" |
    sh
done

mkdir -p data/GSE154778_RAW
mv GSM* data/GSE154778_RAW
