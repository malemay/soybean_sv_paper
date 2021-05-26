#!/bin/bash

ulimit -n 2048

# Getting the path to the executables from the command line
kmc=$1
btools=$2

# DEPENDENCY : breakpoint_refinement_analysis/illumina_ids.txt
for i in $(cat ../illumina_ids.txt)
do
  mkdir $i
  cd $i
 # DEPENDENCY : trimmed Illumina fastq files
  for j in $(ls ../../illumina_data/trimmed_fastq/${i}/${i}*.gz)
  do
    echo $j >> file_list.txt
  done 
 
  mkdir tmpdir
  $kmc -t4 -m12 -sm -k55 -ci1 @file_list.txt $i tmpdir
  $btools makeBloom -p 4 -k $i
  cd ..

done

