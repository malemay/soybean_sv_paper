#!/bin/bash

# Might be necessary to include Python 3 on the path for Nanoplot to run properly
# Loading the python module to run Nanoplot
# module load python/3.5

# Getting the path to the Nanoplot executable from file
NanoPlot=$1

# Iterating over all the 17 lines 
# DEPENDENCY : utilities/line_ids.txt
for i in $(tail -n+2 ../utilities/line_ids.txt | cut -f1)
do
	# DEPENDENCY : Nanopore reads aligned with NGMLR and sorted with Samtools
	$NanoPlot --bam ../nanopore_data/${i}_porechopped_aligned.sort.bam --loglength --plots dot kde hex pauvre --N50 -o $i
done

