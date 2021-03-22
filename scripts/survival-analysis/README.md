# Survival analysis of cancer patients

The codes to generate FigXXX

## Dependencies

- wget
- gzip
- R

## Usage

```bash

git clone https://github.com/bioinfo-tsukuba/SSD.git

# download and preprocess data
./SSD/scripts/survival-analysis/survival.sh

# calculate generalized wilcoxon test
./SSD/scripts/survival-analysis/survival_pval.R

# plot Kaplan-Meier curve
./SSD/scripts/survival-analysis/survival_plot.R

```
