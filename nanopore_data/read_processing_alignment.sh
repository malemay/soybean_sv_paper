#!/bin/bash

# Getting the paths for the various executables from the command line
porechop=$1
ngmlr=$2
samtools=$3

# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# DEPENDENCY : utilities/flowcell_names.txt
# DEPENDENCY : nanopore_data/AC2001_1.fastq.gz
# DEPENDENCY : nanopore_data/AC2001_2.fastq.gz
# DEPENDENCY : nanopore_data/ALTA.fastq.gz
# DEPENDENCY : nanopore_data/MAPLE_ISLE_1.fastq.gz
# DEPENDENCY : nanopore_data/MAPLE_ISLE_2.fastq.gz
# DEPENDENCY : nanopore_data/MAPLE_PRESTO.fastq.gz
# DEPENDENCY : nanopore_data/OAC_09_35C.fastq.gz
# DEPENDENCY : nanopore_data/OAC_CARMAN.fastq.gz
# DEPENDENCY : nanopore_data/OAC_DRAYTON.fastq.gz
# DEPENDENCY : nanopore_data/OAC_EMBRO.fastq.gz
# DEPENDENCY : nanopore_data/OAC_LAKEVIEW_1.fastq.gz
# DEPENDENCY : nanopore_data/OAC_LAKEVIEW_2.fastq.gz
# DEPENDENCY : nanopore_data/OAC_MADOC.fastq.gz
# DEPENDENCY : nanopore_data/OAC_OXFORD.fastq.gz
# DEPENDENCY : nanopore_data/OAC_PETREL_2.fastq.gz
# DEPENDENCY : nanopore_data/OAC_PRUDENCE.fastq.gz
# DEPENDENCY : nanopore_data/OAC_STRATFORD.fastq.gz
# DEPENDENCY : nanopore_data/OT09-03.fastq.gz
# DEPENDENCY : nanopore_data/QS5091.50J.fastq.gz
# DEPENDENCY : nanopore_data/ROLAND.fastq.gz

# Running Porechop using GNU parallel (4 x 4 cores)
cat ../utilities/flowcell_names.txt | parallel -j4 "$porechop --threads 4 --discard_middle -i {}.fastq.gz -o {}_porechopped_sequences.fastq"

# Aligning all reads with NGMLR
cat ../utilities/flowcell_names.txt | parallel -j4 "$ngmlr -r $refgenome -q ${i}_porechopped_sequences.fastq -o ${i}_porechopped_aligned.sam -x ont"

# Converting the sam files into bam
cat ../utilities/flowcell_names.txt | parallel -j4 "$samtools view -b {}_porechopped_aligned.sam > {}_porechopped_aligned.bam"

# Sorting the bam files
cat ../utilities/flowcell_names.txt | parallel -j4 "$samtools sort {}_porechopped_aligned.bam > {}_porechopped_aligned.sort.bam"

# Indexing the sorted bam files
cat ../utilities/flowcell_names.txt | parallel -j4 "$samtools index {}_porechopped_aligned.sort.bam"

# Iterating over the three samples that were sequenced using two different libraries to merge their bam files together
for i in AC2001 MAPLE_ISLE OAC_LAKEVIEW
do
	# Merging the bam files from libraries 1 and 2 into a single sorted file with samtools merge
	$samtools merge ${i}_porechopped_aligned.sort.bam \
		${i}_1_porechopped_aligned.sort.bam \
		${i}_2_porechopped_aligned.sort.bam

	# Indexing the newly created file
	$samtools index ${i}_porechopped_aligned.sort.bam
done

# DEPENDENCY : nanopore_data/AC2001_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/ALTA_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/MAPLE_ISLE_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/MAPLE_PRESTO_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_09_35C_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_CARMAN_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_DRAYTON_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_EMBRO_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_LAKEVIEW_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_MADOC_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_OXFORD_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_PETREL_2_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_PRUDENCE_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_STRATFORD_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OT09-03_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/QS5091.50J_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/ROLAND_porechopped_aligned.sort.bam

