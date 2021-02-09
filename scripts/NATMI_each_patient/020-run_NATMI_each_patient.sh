#!/bin/bash

# Change working directory
cd NATMI

# Definition of input and output paths
indir="../processed_data/NATMI_each_patient"
outdir_base="../results/NATMI_each_patient/ExtractEdges"

# Make output base directory
mkdir -p $outdir_base

# Get patient sample ids
samples=`ls ${indir}/*_annot.csv | awk '{gsub(/.+\//,""); gsub("_annot.csv", ""); print}'`

# For each patient. run ExtractEdges.py
for sample in $samples
do
    echo $sample
    emFile="${sample}_expression.csv"
    annFile="${sample}_annot.csv"
    out="${outdir_base}/${sample}"

    python ExtractEdges.py --species human --emFile ${emFile} --annFile ${annFile} --interDB lrc2p --out ${out}

done

