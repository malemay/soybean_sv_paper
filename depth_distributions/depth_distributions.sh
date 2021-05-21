#!/bin/bash

# Getting the path to the samtools executable from the command line
samtools=$1

# Running the bash/awk script on it
# DEPENDENCY : utilities/line_ids.txt
# DEPENDENCY : scripts/depth_distribution.sh
# DEPENDENCY : Oxford Nanopore reads aligned with Sniffles
tail -n+2 ../utilities/line_ids.txt | cut -f1 | \
	parallel -j6 "../scripts/depth_distribution.sh ../nanopore_data/{}_porechopped_aligned.sort.bam {}_depth.txt $samtools"

