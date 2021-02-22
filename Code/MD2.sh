#MD2 run NATMI_ExtractEdges

#!/bin/bash


# Change working directory
cd NATMI

# Definition of input and output paths
indir="../processed_data/MD2"
outdir_base="../results/MD2/ExtractEdges"

# Make output base directory
mkdir -p $outdir_base

# Get patient sample ids
samples=`ls ${indir}/*_annot.csv | awk '{gsub(/.+\//,""); gsub("_annot.csv", ""); print}'`

# For each patient. run ExtractEdges.py

for sample in $samples
do
    echo $sample
    emFile="$indir/${sample}_expression.csv"
    annFile="$indir/${sample}_annot.csv"
    out="${outdir_base}/${sample}"

    /Users/sayakasuzuki/opt/anaconda3/envs/ssd/bin/python ExtractEdges.py --species human --emFile ${emFile} --annFile ${annFile} --interDB lrc2p --out ${out}

done

#run NATMI_DiffEdges.sh

/Users/sayakasuzuki/opt/anaconda3/envs/ssd/bin/python DiffEdges.py --refFolder "../results/MD2/ExtractEdges/P_low_grade" --targetFolder "../results/MD2/ExtractEdges/P_high_grade" --interDB lrc2p --out "../results/MD2/DiffEdges"

#run NATMI_VisInteractions.py
#Visualise cell-connectivity-summary networks from the results
#Visualize the low grade
/Users/sayakasuzuki/opt/anaconda3/envs/ssd/bin/python VisInteractions.py --sourceFolder "../results/MD2/ExtractEdges/P_low_grade" --interDB lrc2p --weightType mean --detectionThreshold 0.2 --plotFormat pdf --drawNetwork y --plotWidth 12 --plotHeight 10 --layout kk --fontSize 8 --edgeWidth 0 --maxClusterSize 0 --clusterDistance 1

#Visualize the high grade
/Users/sayakasuzuki/opt/anaconda3/envs/ssd/bin/python VisInteractions.py --sourceFolder "../results/MD2/ExtractEdges/P_high_grade" --interDB lrc2p --weightType mean --detectionThreshold 0.2 --plotFormat pdf --drawNetwork y --plotWidth 12 --plotHeight 10 --layout kk --fontSize 8 --edgeWidth 0 --maxClusterSize 0 --clusterDistance 1

#Visualize the DiffEdges
/Users/sayakasuzuki/opt/anaconda3/envs/ssd/bin/python VisInteractions.py --sourceFolder "../results/MD2/DiffEdges" --interDB lrc2p --weightType mean --detectionThreshold 0.2 --plotFormat pdf --drawNetwork y --plotWidth 12 --plotHeight 10 --layout kk --fontSize 8 --edgeWidth 0 --maxClusterSize 0 --clusterDistance 1

#draw cluster pair 
/Users/sayakasuzuki/opt/anaconda3/envs/ssd/bin/python VisInteractions.py --sourceFolder "../results/MD2/DiffEdges" --interDB lrc2p --drawClusterPair y
