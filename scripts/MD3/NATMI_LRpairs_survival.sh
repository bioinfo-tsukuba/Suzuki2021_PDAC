#動作確認
Rscript ./scripts/MD3/lrpair_pval.R


#NATMI LR pair survival analysis

#PAAD-US patients
Rscript ./scripts/MD3/lrpair_pval.R \
  data/ICGC/survival_PAAD-US.csv.gz \
  data/NATMI_lrpair.csv > results/NATMI_lrpair_result.csv


#PACA-AU patients
Rscript ./scripts/MD3/lrpair_pval.R \
  data/ICGC/survival_PACA-AU.csv.gz \
  data/NATMI_lrpair.csv > results/NATMI_lrpair_result_AU.csv

  #PACA-CA patients
Rscript ./scripts/MD3/lrpair_pval.R \
  data/ICGC/survival_PACA-CA.csv.gz \
  data/NATMI_lrpair.csv > results/NATMI_lrpair_result_CA.csv
