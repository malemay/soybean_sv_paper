#!/bin/bash

# Getting the path to the FLASH executable from the command line
flash=$1

# Moving the the flash directory
cd flash_merging

# Processing the samples that have a read length of 100
for i in $(seq 1001 1023 | xargs -I {} echo CAD{})
do
  mkdir -p $i
  fastq_dir=../../../illumina_data/trimmed_fastq/${i}
  # DEPENDENCY : paired reads output by bbduk
  $flash -z -M 100 -d ${i} -t 1 ${fastq_dir}/${i}_R1_trimmed.fastq.gz ${fastq_dir}/${i}_R2_trimmed.fastq.gz
done

# Processing the samples that have a read length of 101
for i in $(seq 1025 1087 | xargs -I {} echo CAD{})
do
  mkdir -p $i
  fastq_dir=../../../illumina_data/trimmed_fastq/${i}
  # DEPENDENCY : paired reads output by bbduk
  $flash -z -M 101 -d ${i} -t 1 ${fastq_dir}/${i}_R1_trimmed.fastq.gz ${fastq_dir}/${i}_R2_trimmed.fastq.gz
done

# Processing the samples that have a read length of 125
for i in $(seq 1088 1095 | xargs -I {} echo CAD{})
do
  mkdir -p $i
  fastq_dir=../../../illumina_data/trimmed_fastq/${i}
  # DEPENDENCY : paired reads output by bbduk
  $flash -z -M 101 -d ${i} -t 1 ${fastq_dir}/${i}_R1_trimmed.fastq.gz ${fastq_dir}/${i}_R2_trimmed.fastq.gz
done

# Processing the samples that have a read length of 126
for i in $(seq 1096 1103 | xargs -I {} echo CAD{})
do
  mkdir -p $i
  fastq_dir=../../../illumina_data/trimmed_fastq/${i}
  # DEPENDENCY : paired reads output by bbduk
  $flash -z -M 101 -d ${i} -t 1 ${fastq_dir}/${i}_R1_trimmed.fastq.gz ${fastq_dir}/${i}_R2_trimmed.fastq.gz
done

