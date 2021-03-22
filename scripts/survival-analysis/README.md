# Survival analysis of cancer patients

The codes to generate FigXXX

## Dependencies

- wget
- gzip
- R

## Usage

### 1. Download and preprocess data

```
./SSD/scripts/survival-analysis/survival.sh
```
Output file is `data/ICGC/survival_*.csv.gz`.

| id       | sex  | status | time | gene  | exp     |
| -------- | ---- | ------ | ---- | ----- | ------- |
| DO221539 | male | alive  | 1733 | A2M   | 1063.3  |
| DO221539 | male | alive  | 1733 | AANAT | 0       |
| DO221539 | male | alive  | 1733 | ABCA1 | 144.511 |

### 2. Calculate generalized wilcoxon test

```
./SSD/scripts/survival-analysis/survival_pval.R [survival data] [LR pair]
```

### 3. Plot Kaplan-Meier curve

```
./SSD/scripts/survival-analysis/survival_plot.R  [survival data] [LR pair]
```
