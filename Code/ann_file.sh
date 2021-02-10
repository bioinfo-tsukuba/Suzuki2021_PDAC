# input を output に変換したい
# なぜなら こういう形式にしたいから https://github.com/asrhou/NATMI/blob/master/toy.sc.ann.txt
input=SSD/data/WeiLin_pdac10/PDAC10_meta.csv
output=SSD/processed_data/PDAC10_metaFile.tsv
# どんな入力ファイルか確認しよう！
head $input
# 変換するよ！
awk -F"," 'BEGIN{OFS="\t"; print "barcode","annotation"} NR>1{print $1,$3}' $input > $output
# どんな出力ファイルになったか確認しよう
#（ https://github.com/asrhou/NATMI/blob/master/toy.sc.ann.txt っぽくなったかな？？ ）
head $output