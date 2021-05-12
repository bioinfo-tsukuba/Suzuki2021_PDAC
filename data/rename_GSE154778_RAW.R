#################################
# このスクリプトは、data/GSE154778_RAWの
# Seurat パッケージの `Read10X()` 関数が読み込めるようする。
# そのために、サンプルごとにフォルダを作成し、
# その中にbarcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gzがある状態にする。
# また、元のファイル名と変更されたファイル名の対応表も作成する。
#################################

# ライブラリを読み込む
library(tidyverse)

# 出力ディレクトリを作成する
outdir <- "data/GSE154778_RAW_renamed"
dir.create(outdir, recursive = TRUE)

# `data/GSE154778_RAW` に含まれるファイルのパスを取得する
filenames <- list.files("data/GSE154778_RAW", "*.gz", full.names = TRUE)

# 元のファイル名が含まれたデータフレームを作成する
df1 <- tibble(filename_original = filenames)

# ファイル名からサンプル名をパースし、新しいファイル名を filename_renamed 列に入れる
df1 %>%
  mutate(Sample_ID = gsub("_[^_]+\\.[^_]+\\.gz$", "", basename(filename_original))) %>%
  mutate(filename_renamed = file.path(outdir, Sample_ID, gsub("^.+_", "", basename(filename_original)))) %>%
  mutate(filename_renamed = gsub("genes.tsv.gz", "features.tsv.gz", filename_renamed)) ->
df1

# df1 で作成された「元のファイル名と変更されたファイル名の対応表」に従い、
# ファイルを名前を変えつつコピーする
for (i in 1:nrow(df1)) {
  newdir <- dirname(df1[i, ]$filename_renamed)
  if (!dir.exists(newdir)) {
    dir.create(newdir, recursive = TRUE)
  }
  file.copy(df1[i, ]$filename_original, df1[i, ]$filename_renamed)
}

# 「元のファイル名と変更されたファイル名の対応表」を保存する
df1 %>%
  write_tsv(file.path(outdir, "filenames.tsv"))

# SessionInfo()
sessionInfo()
