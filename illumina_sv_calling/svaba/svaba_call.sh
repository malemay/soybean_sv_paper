#!/bin/bash

# Getting the path to the SvABA executable from the command line
svaba=$1

# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# DEPENDENCY : utilities/all_lines.txt
# DEPENDENCY : Illumina reads aligned with BWA
# Running 4 svaba tasks of 4 threads each in parallel
cat ../../utilities/all_lines.txt | parallel -j 4 "$svaba run -t ../../illumina_data/aligned_reads/{}/{}_all.sort.bam -p 4 -a {} -G $refgenome --germline -I -L 6 --verbose 4"

