# Survival analysis of cancer patients

The codes to generate FigXXX

## Dependencies

- wget
- gzip
- R

## Usage

### Download and preprocess data

```sh
./scripts/Fig1/00-format_ICGC.sh
```

The script downloads and formats ICGC PDAC cohort data.  

Output files are followings.
- data/ICGC/survival_PAAD-US.csv.gz
- data/ICGC/survival_PACA-AU.csv.gz
- data/ICGC/survival_PACA-CA.csv.gz

The formats of `data/ICGC/survival_*.csv.gz`:point_down:

| id       | sex  | status | time | gene  | exp     |
| -------- | ---- | ------ | ---- | ----- | ------- |
| DO221539 | male | alive  | 1733 | A2M   | 1063.3  |
| DO221539 | male | alive  | 1733 | AANAT | 0       |
| DO221539 | male | alive  | 1733 | ABCA1 | 144.511 |

The `exp` represents CPM normalized gene expression.

### Calculate P-value and HR for each ligand-receptor pair

```sh
./scripts/Fig1/02-pval_ICGC.sh
```

Output file is:
- results/Fig1/LR_Pval_HR.csv

## Extract LR pairs with adj-Pval < 0.1 and consistent HR values (all cohorts shows HR>1 or HR<1)

```sh
Rscript --slave --vanilla ./scripts/Fig1/03-pval_metap.R
```

Output file is:
- results/Fig1/LR_adjPval_meanHR_screened.csv

The file contains 199 LR pairs (67 LRs in HR<1 and 133 LRs in HR>1)
