#!/bin/sh

#==============================================================================
# Initialization
#==============================================================================

set -u
umask 0022
export LC_ALL=C
[ -z "$ZSH_VERSION" ] || setopt shwordsplit interactivecomments

# --- cURL or Wget ----------------------------------------
if type curl >/dev/null 2>&1; then
  CMD_GET="curl -Lso -"
elif type wget >/dev/null 2>&1; then
  CMD_GET="wget --no-check-certificate -qO -"
else
  { echo 'No HTTP-GET/POST command found.'; exit 1; }
fi

# --- gzip ------------------------------------------------
type gzip >/dev/null 2>&1 || { echo 'gzip: command not found.'; exit 1; }

#==============================================================================
# Retrieve gene lists related to ligand-receptor interaction from NATMI
#==============================================================================

mkdir -p data/HGNC

$CMD_GET https://asrhou.github.io/NATMI/ |
  awk 'BEGIN{RS=""} {$1=$1}1' |
  grep '<td class="col1">' |
  sed "s/.*<td class=\"col1\"> //" |
  sed "s/ .*<td class=\"col4\">//" |
  cut -d " " -f 1,2 |
  tr " " "\n" |
  sort -u > data/HGNC/natmi_genes.txt

#==============================================================================
# Retrieve HGNC gene list to convert Ensemble to Gene Symbol
#==============================================================================

mkdir -p data/HGNC

$CMD_GET ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt |
  sed 1d |
  cut -f 2,20 |
  awk 'NF==2' |
  sort |
  join - data/HGNC/natmi_genes.txt |
  awk '{print $2,$1}' |
  sort -t " " > data/HGNC/natmi_ensemble_symbol.txt

cut -d " " -f 1 data/HGNC/natmi_ensemble_symbol.txt | sort | awk '{print $1,"LR"}' > data/HGNC/natmi_ensemble.txt
cut -d " " -f 2 data/HGNC/natmi_ensemble_symbol.txt | sort | awk '{print $1,"LR"}' > data/HGNC/natmi_symbol.txt

#==============================================================================
# Retrieve patients data from ICGC
#==============================================================================

mkdir -p data/ICGC
echo "*" > data/ICGC/.gitignore

for project in PAAD-US PACA-AU PACA-CA; do
cat << FILE |
https://dcc.icgc.org/api/v1/download?fn=/current/Projects/$project/donor.$project.tsv.gz
https://dcc.icgc.org/api/v1/download?fn=/current/Projects/$project/exp_seq.$project.tsv.gz
FILE
  while read -r file; do
    $CMD_GET "$file" > "data/ICGC/${file##*/}"
  done
done

# The number of patients ----------------------------------

find data/ICGC/ -type f |
  grep donor |
  while read -r file; do
    echo "$file"
    gzip -dc "$file" |
      grep -v "icgc_donor_id" |
      wc -l
  done

# CPM normalization by patiants ---------------------------

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

#==============================================================================
# Generate CSV files including "Patient ID", "Sex", "Status", "Survival time", "Gene", "Expression"
#==============================================================================

mkdir -p data/ICGC/

gzip -dc data/ICGC/donor.* |
  grep -v icgc_donor_id |
  cut -f 1,5,6,17,18 |
  tr "\t" "@" |
  awk -F "@" 'BEGIN{OFS=" "} $3=="alive" {$4=$5} {$1=$1}1' |
  cut -d " " -f 1-4 |
  awk 'NF==4' |
  sort -t " " > tmp_donor

cat data/ICGC/exp_seq_* |
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
  cat > data/ICGC/survival.csv

rm tmp_donor
