#!/bin/bash

module load deepTools/3.5.5-foss-2022b
cd /g/krebs/barzaghi/analyses/31.01.23_GenVar_figures/TF-chromatin/revisions_02-06-2025

plotHeatmap \
  -m deepTools_matrix.gz \
  -o deepTools_heatmap.pdf \
  --sortUsingSamples 1 \
  -x "rel distance (bp)" \
  --startLabel "-2500" \
  --endLabel "+2500" \
  --regionsLabel "enhancers" "promoters" "intergenic sites" "CTCF sites" \
  --samplesLabel "DNase" "H3K27ac" "H3K4me1" "H3K4me3" "CTCF" \
  --zMin 0 0 0 0 0 \
  --zMax 40 40 20 80 20