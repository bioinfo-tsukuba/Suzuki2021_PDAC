# Survival analysis of cancer patients

The codes to generate FigXXX

## Dependencies

- wget
- gzip
- R

## Usage

### Download and preprocess data

```bash
./SSD/library/survival.sh
```

ファイルの保存先はこんな感じです:point_right: `data/ICGC/survival_*.csv.gz`.  

`data/ICGC/survival_*.csv.gz`の形式はこんな感じです:point_down:

| id       | sex  | status | time | gene  | exp     |
| -------- | ---- | ------ | ---- | ----- | ------- |
| DO221539 | male | alive  | 1733 | A2M   | 1063.3  |
| DO221539 | male | alive  | 1733 | AANAT | 0       |
| DO221539 | male | alive  | 1733 | ABCA1 | 144.511 |


### 各細胞ごとに生存時間解析を行うコマンド

#### P-value

```bash
Rscript library/CCI_pval.R \
  data/ICGC/survival_PAAD-US.csv.gz \
  tests/MD3/data/CCI.csv > results_CCI.csv
```

入力ファイル`tests/MD3/data/CCI.csv`の形式はこんな感じです:point_down:  
| CCI      | LR          |
| -------- | ----------- |
| DC->EMT  | VIM->CD44   |
| CAF->EMT | TIMP1->CD63 |
| TIL->EMT | VIM->CD44   |


出力ファイル`results_CCI.csv`の形式はこんな感じです:point_down:  
| CCI       | pval | median_day_low | median_day_high | diff_day |
| --------- | ---- | -------------- | --------------- | -------- |
| Endo->EMT | 0    | 180            | 188             | 8        |
| CAF->EMT  | 0    | 180            | 188             | 8        |
| TAM->EMT  | 0    | 180            | 188             | 8        |

#### Plot

```bash
Rscript ./library/CCI_plot.R \
  data/ICGC/survival_PAAD-US.csv.gz \
  tests/MD3/data/CCI.csv
```

### 各リガンドレセプターごとに生存時間解析を行うコマンド

#### P-value

```bash
Rscript ./library/LR_pval.R \
  data/ICGC/survival_PAAD-US.csv.gz \
  tests/MD3/data/LR.csv > results_LR.csv
```

入力ファイル`tests/MD3/data/LR.csv`の形式はこんな感じです:point_down:  

| LR           |
| ------------ |
| VIM->CD44    |
| TIMP1->CD63  |
| COL1A2->CD44 |

出力ファイル`results_LR.csv`の形式はこんな感じです:point_down:  

| LR            | pval                  | median_day_low | median_day_high | diff_day |
| ------------- | --------------------- | -------------- | --------------- | -------- |
| ICOSLG->CTLA4 | 3.834185013928959e-8  | 165.5          | 216             | 50.5     |
| ICAM4->ITGAL  | 4.794087982151751e-8  | 185.5          | 181             | -4.5     |
| ICOSLG->CD28  | 1.1746727190953266e-7 | 181            | 216             | 35       |

#### Plot

```bash
Rscript ./library/LR_plot.R \
  data/ICGC/survival_PAAD-US.csv.gz \
  tests/MD3/data/LR_top10.csv
```