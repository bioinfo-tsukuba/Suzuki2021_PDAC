# Survival meta analysis

## やったこと

`run.sh` の実行

1. Run survival analysis: "PAAD-US", "PACA-AU", "PACA-CA" の３コホートのデータに対し、`library/LR_pval.R` を実行する。
2. Run metap analysis: 1の結果に対し、`metap.R` を実行する。
   1. Combine p-values using the sum of p (Edgington's) method

## 結果の場所

`results/`

## `metap.R` の出力 (`results_LR_meta.csv`) の解析方法

```r
# BH法による q-value < 0.1 && ３つの cohort での HRが < 1 である行（LRペアを抽出する）
df_LR_meta %>% 
  filter(meta_qval_BH < 0.1, `hr_PAAD-US` < 1, `hr_PACA-AU` < 1, `hr_PACA-CA` < 1)

# Holm法による q-value < 0.1 && ３つの cohort での HRが < 1 である行（LRペアを抽出する）
df_LR_meta %>% 
  filter(meta_qval_holm < 0.1, `hr_PAAD-US` < 1, `hr_PACA-AU` < 1, `hr_PACA-CA` < 1)

# Storey法による q-value < 0.1 && ３つの cohort での HRが < 1 である行（LRペアを抽出する）
df_LR_meta %>% 
  filter(meta_qval_storey < 0.1, `hr_PAAD-US` < 1, `hr_PACA-AU` < 1, `hr_PACA-CA` < 1)
```

## 更新履歴

- 2021/04/17  ver 0.1.0
