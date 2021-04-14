# Run NATMI
cd NATMI 

conda activate ssd

path_emFile=../SSD/processed_data/expression_cluster1.csv
out=2021_02_09
path_annFile=../SSD/processed_data/PDAC10_metaFile.tsv
python ExtractEdges.py --species human --emFile ${path_emFile} --annFile ${path_annFile} --interDB lrc2p --coreNum 4 --out $out

python VisInteractions.py --sourceFolder  $out --interDB lrc2p --weightType mean --detectionThreshold 0.2 --drawNetwork y --plotWidth 4 --plotHeight 4 --layout circle --fontSize 15 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 0.6

python VisInteractions.py --sourceFolder $out --drawClusterPair y --keepTopEdge 15


# to visialize the cell-to-cell communication network via CD2-CD48 pairs
python VisInteractions.py --sourceFolder $out --drawLRNetwork CD2 CD48 --plotWidth 4 --plotHeight 4 --layout circle --fontSize 15 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 0.6
