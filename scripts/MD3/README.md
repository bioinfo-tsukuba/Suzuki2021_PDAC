# Survival analysis of cancer patients

The codes to generate FigXXX

## Dependencies

- wget
- gzip
- R

## Usage

### 1. Download and preprocess data

```
./SSD/scripts/MD3/survival.sh
```

ファイルの保存先はこちらです→ `data/ICGC/survival_*.csv.gz`.  
`data/ICGC/survival_*.csv.gz`の形式↓

| id       | sex  | status | time | gene  | exp     |
| -------- | ---- | ------ | ---- | ----- | ------- |
| DO221539 | male | alive  | 1733 | A2M   | 1063.3  |
| DO221539 | male | alive  | 1733 | AANAT | 0       |
| DO221539 | male | alive  | 1733 | ABCA1 | 144.511 |


### 各細胞ごとに生存時間解析を行うコマンド

```
./SSD/scripts/MD3/celltype_lrpair_pval.R data/ICGC/survival_PAAD-US.csv.gz celltype_and_lrpair.csv > celltype_lrpair_result.csv
```

入力ファイル`celltype_and_lrpair.csv`の形式↓  
| celltype_pair | ligandreceptor_pair |
| ------------- | ------------------- |
| DC->EMT       | VIM->CD44           |
| CAF->EMT      | TIMP1->CD63         |
| TIL->EMT      | VIM->CD44           |

出力ファイル`celltype_lrpair_result.csv`の形式↓  
| celltype_pair | pval | median_day_low | median_day_high | diff_day |
| ------------- | ---- | -------------- | --------------- | -------- |
| Endo->EMT     | 0    | 180            | 188             | 8        |
| CAF->EMT      | 0    | 180            | 188             | 8        |
| TAM->EMT      | 0    | 180            | 188             | 8        |


### 各リガンドレセプターごとに生存時間解析を行うコマンド

```
./SSD/scripts/MD3/lrpair_pval.R data/ICGC/survival_PAAD-US.csv.gz lrpair.csv > lrpair_result.csv
```

入力ファイル`lrpair.csv`の形式↓  

| ligandreceptor_pair |
| ------------------- |
| VIM->CD44           |
| TIMP1->CD63         |
| COL1A2->CD44        |

出力ファイル`lrpair_result.csv`の形式↓  

| lr_pair       | pval                  | median_day_low | median_day_high | diff_day |
| ------------- | --------------------- | -------------- | --------------- | -------- |
| ICOSLG->CTLA4 | 3.834185013928959e-8  | 165.5          | 216             | 50.5     |
| ICAM4->ITGAL  | 4.794087982151751e-8  | 185.5          | 181             | -4.5     |
| ICOSLG->CD28  | 1.1746727190953266e-7 | 181            | 216             | 35       |
