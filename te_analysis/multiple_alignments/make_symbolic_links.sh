#!/bin/bash

# Creating symbolic links to the bam files locally to facilitate the analysis
# DEPENDENCY : utilities/line_ids.txt
for i in $(tail -n+2 ../../utilities/line_ids.txt | cut -f1) 
do
	ln -s ../../nanopore_data/${i}_porechopped_aligned.sort.bam
	ln -s ../../nanopore_data/${i}_porechopped_aligned.sort.bam.bai
done

