# Survival meta analysis

## やったこと

- LRペアの生存時間解析の３コホートの結果のメタ解析の結果の全体像を示す図を作成した。
  - 散布図
    - 点：LRペア
    - x軸：HRの平均（算術平均、幾何平均）
    - y軸：メタ解析のp-valueのStorey法の-log10(adjusted p-value)

### 実行方法

```bash
bash scripts/survival_meta_analysis_plot/run.sh
```


## 更新履歴

- 2021/05/11  ver 1.0.0
