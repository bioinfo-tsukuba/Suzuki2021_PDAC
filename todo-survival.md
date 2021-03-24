# 生存時間解析

## ToDO


## DONE

### 2021-03-23

+ [x] 出力の`data/ICGC/survival_*.csv.gz`の並び順が異なっていたので, 目的のもの`id,sex,status,time,gene,exp`に変更しました.

### 2021-03-22

+ [x] ICGCからPrimaryのみにする
+ [x] ICGCをコホートごとにする

### 2021-03-19
+ [x] LRペア全パターンの生存時間
  + => `results/MD3/survival_all_LR.csv`に保存しました

+ [x] リガンドレセプターに関わる遺伝子群を取ってくる
+ [x] ICGCから膵臓がんのデータを取ってくる
+ [x] 患者さんID, 性別, 生死, 観測期間, 遺伝子名, 遺伝子発現の6コラムのCSVファイルをつくる
+ [x] 遺伝子名を統一する（ENSEMBLEをGene Symbolへ）
+ [x] 発現データからリガンドレセプターに関わる遺伝子のみを抽出する
+ [x] (哀しみ) プロジェクトごとに正規化リードカウントが全然違う値なので, 生のリードカウントをCPM正規化する
+ [x] survminerパッケージを用いて生存曲線をかく
+ [x]  generalized Wilcoxon testを出力する
+ [x]  category(Appeared, Disappered, Up, Down)ごとのLRペアについてP-valを出力する
+ [x] delta_edge_specificity_weight 以外の指標に変えると結果がどうなるか？
  + -> **expression_weightでは"CD44"が関わっていそうで面白いです.生存曲線でCD44はどうなっているのが興味があります.**
+ [x] ICGCはいろんなステージが混ざっているので、Primary vs. Metastasis でやるとまた異なる結果になる？
  + -> **すべてPrimaryまたはUnknownのようです**
