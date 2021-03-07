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
https://dcc.icgc.org/api/v1/download?fn=/current/Projects/$project/specimen.$project.tsv.gz
FILE
  while read -r file; do
    echo "$file"
    [ -s "data/ICGC/${file##*/}" ] || $CMD_GET "$file" > "data/ICGC/${file##*/}"
    # in case of failed first try...
    [ -s "data/ICGC/${file##*/}" ] || $CMD_GET "$file" > "data/ICGC/${file##*/}"
  done
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
# Summarize specimen info
#==============================================================================

cat data/ICGC/nationwidechildrens.org_clinical_patient_paad.txt |
  grep -v -e bcr_patient_barcode -e CDE_ID |
  tr " ,\t" "@@ " |
  cut -d " " -f2,24,32,70 |
  sed "s/Stage@//" |
  sed "s/\[Not@Available\]/NA/" |
  sed "s/\[Discrepancy\]/Discrepancy/" |
  sort -t " " |
  cat > tmp_tcga

gzip -dc data/ICGC/specimen.PAAD-US.tsv.gz |
  grep -v "icgc_donor_id" |
  tr " " "_" |
  awk '{print $5,$2,$4}' |
  sort -t " " |
  join - tmp_tcga |
  # project, donor id, tumor, grade, stage
  awk '{print $2,$3,$6,$4,$5}' |
  cat > tmp_tcga_patientinfo

gzip -dc data/ICGC/specimen.PACA-CA.tsv.gz data/ICGC/specimen.PACA-AU.tsv.gz |
  grep -v "icgc_donor_id" |
  tr " " "_" |
  # project, donor id, tumor, grade, stage
  cut -f 2,5,20,22,25 |
  awk 'NF == 5' |
  cat - tmp_tcga_patientinfo |
  sort -u |
  # Select PDAC
  grep -e "Pancreatic_Ductal_Adenocarcinoma" -e "8140/3" -e "8500/3" -e "8021/3" -e "8035/3" |
  # Format Grade
  awk '
    $4=="1_-_Well_differentiated" {$4="G1"}
    $4=="2_-_Moderately_differentiated" {$4="G2"}
    $4=="3_-_Poorly_differentiated" {$4="G3"}
    $4=="4_-_Undifferentiated" {$4="G4"}
    $4=="Moderate_to_Poor" {$4="G2-G3"}
    $4=="Moderate_to_Poorly_differentiated" {$4="G2-G3"}
    $4=="Moderately_differentiated" {$4="G2"}
    $4=="Not_specified/Unknown" {$4="X"}
    $4=="Other_-_Status_Post_Therapy" {$4="X"}
    $4=="Poorly_differentiated" {$4="G3"}
    $4=="Poorly_differentiated_to_Undifferentiated" {$4="G3-G4"}
    $4=="Undifferentiated" {$4="G4"}
    $4=="Well_differentiated" {$4="G1"}
    $4~"X" {$4="X"}1
  ' |
  # Annotate Primary or Metastasis
  awk '$NF ~ "M1" || $NF ~ "IV" {$(NF+1)="Metastasis"} {$(NF+1)="Primary"}1' |
  awk '{print $2,$1,$4,$6,$5}' |
  cat > tmp_patients_info

rm tmp_tcga tmp_tcga_patientinfo

#==============================================================================
# Generate CSV files including "Patient ID", "Sex", "Status", "Survival time", "Gene", "Expression"
#==============================================================================

mkdir -p data/ICGC/

cut -d "," -f 1 tmp_patients_info | sort -u > tmp_patients_id

gzip -dc data/ICGC/donor.* |
  grep -v icgc_donor_id |
  cut -f 1,5,6,17,18 |
  tr "\t" "@" |
  awk -F "@" 'BEGIN{OFS=" "} $3=="alive" {$4=$5} {$1=$1}1' |
  cut -d " " -f 1-4 |
  awk 'NF==4' |
  sort -t " " |
  join - tmp_patients_id > tmp_donor

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

rm tmp_patients_id tmp_donor

#==============================================================================
# Output patients info
#==============================================================================

cut -d, -f 1 data/ICGC/survival.csv | sort -u > tmp_donor_w_rnaseq

sort -t " " tmp_patients_info |
  join - tmp_donor_w_rnaseq |
  awk 'BEGIN{OFS=","; print "ID", "Project", "Grade", "Primary_Metastasis", "Stage"}1' |
  tr " " "," |
  cat > results/MD3/patients_info.csv

rm tmp_patients_info tmp_donor_w_rnaseq
