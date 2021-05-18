#!/bin/bash

bwa=$1
samtools=$2
bamaddrg=$3

# Moving to the alignment directory
cd aligned_reads

# DEPENDENCY : utilities/all_lines.txt
for i in $(cat ../../utilities/all_lines.txt)
do
	mkdir $i
	cd $i
	fastq_path=../../trimmed_fastq/${i}

        # DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
	refgenome=../../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

        # Aligning the paired reads
	$bwa mem -t 8 $refgenome ${fastq_path}/${i}_R1_trimmed.fastq.gz ${fastq_path}/${i}_R2_trimmed.fastq.gz | $samtools view -bSh - > ${i}_paired.bam

        # Aligning the reads that were unpaired after trimming
	$bwa mem -t 8 $refgenome ${fastq_path}/${i}_sing_trimmed.fastq.gz | $samtools view -bSh - > ${i}_unpaired.bam

        # Sorting the reads
	$samtools sort -@ 8 -o ${i}_paired.sort.bam  ${i}_paired.bam
	$samtools sort -@ 8 -o ${i}_unpaired.sort.bam  ${i}_unpaired.bam

        # Merging the paired and unpaired reads
	$samtools merge ${i}_tmp.sort.bam ${i}_paired.sort.bam ${i}_unpaired.sort.bam

        # Indexing the reads
	$samtools index ${i}_tmp.sort.bam 

	# Removing the files that are no longer needed
	rm ${i}_paired.bam ${i}_unpaired.bam ${i}_paired.sort.bam ${i}_unpaired.sort.bam

        # Back to the main directory
	cd ..
done

# Adding read group for each sample
cat ../../utilities/all_lines.txt | parallel -j 12 "$bamaddrg -b {}/{}_tpm.sort.bam -s {} > {}/{}_all.sort.bam"

# Indexing the samples
cat ../../utilities/all_lines.txt | parallel -j 12 "$samtools index {}/{}_all.sort.bam"

# Removing the temporary files
rm */*tmp.sort.bam*

