# 生存時間解析

+ [x] リガンドレセプターに関わる遺伝子群を取ってくる
+ [x] ICGCから膵臓がんのデータを取ってくる
+ [x] 患者さんID, 性別, 生死, 観測期間, 遺伝子名, 遺伝子発現の6コラムのCSVファイルをつくる
+ [x] 遺伝子名を統一する（ENSEMBLEをGene Symbolへ）
+ [x] 発現データからリガンドレセプターに関わる遺伝子のみを抽出する
+ [x] (哀しみ) プロジェクトごとに正規化リードカウントが全然違う値なので, 生のリードカウントをCPM正規化する
+ [x] survminerパッケージを用いて生存曲線をかく
+ [x]  generalized Wilcoxon testを出力する
+ [x]  category(Appeared, Disappered, Up, Down)ごとのLRペアについてP-valを出力する

+ [ ] 生存時間曲線での「LRペアの発現量」の計算にzスコアなどをかませると結果が変わるか？
+ [ ] delta_edge_specificity_weight 以外の指標に変えると結果がどうなるか？

+ [x] ICGCはいろんなステージが混ざっているので、Primary vs. Metastasis でやるとまた異なる結果になる？
  + -> **すべてPrimaryまたはUnknownのようです**
```sh
gzip -dc data/ICGC/donor.*.tsv.gz | tr "\t" "," | cut -d "," -f 13 | sort | uniq -c
```

+ [ ] https://cibersortx.stanford.edu/ のような細胞型ごとの発現量を推定して生存時間解析をやると異なる結果になる？

## 余裕があれば

+ [ ] ポジコンとなる遺伝子で生存曲線を評価する
+ [ ] `survival.R`の引数にsurvival.csvとgenelist.txtの2つを加える
