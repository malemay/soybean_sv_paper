#!/bin/bash

# Getting the paths to the executables from the command line
Rscript=$1
plink=$2
phylip=$3

# Using plink to generate distance matrices from SVs
# DEPENDENCY : structure_analysis/sv_pca/SV_PCA (because of the files in PLINK format created by the PCA script)
$plink --bfile ../sv_pca/sv_pca_input --distance square 1-ibs --out sv_distance

# Converting the distance matrix to a format that PHYLIP can use
$Rscript -e "source('convert_to_phylip.R') ; to_phylip('sv_distance.mdist')"

# Generating the input commands for PHYLIP and running the neighbour-joining tree analysis
echo "sv_distance.mdist.phyl" > sv_phylip_params.txt
echo "Y" >> sv_phylip_params.txt

$phylip < sv_phylip_params.txt > sv_phylip.log

# Changing the name of the output tree
mv outtree sv_outtree
rm outfile

