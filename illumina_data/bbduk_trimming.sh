#!/bin/bash

# Getting the path to the bbduk executable from the command line
bbduk=$1

# Moving to the trimming directory
cd trimmed_fastq

# DEPENDENCY : utilities/all_lines.txt
for i in $(cat ../../utilities/all_lines.txt)
do
  mkdir -p $i
  cd $i
  # DEPENDENCY : Raw Illumina paired-end reads
  $bbduk in1=../../raw_fastq/${i}_R1.fastq.gz in2=../../raw_fastq/${i}_R2.fastq.gz \
      out1=${i}_R1_trimmed.fastq.gz out2=${i}_R2_trimmed.fastq.gz outs=${i}_sing_trimmed.fastq.gz \
      literal=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      stats=${i}_stats_trimmed.txt \
      ktrim=r k=23 mink=11 hdist=1 tpe tbo \
      threads=1 \
      qtrim=r trimq=10 minlength=35 minavgquality=15
  cd ..
done

