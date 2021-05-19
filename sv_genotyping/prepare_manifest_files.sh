#!/bin/bash

# Getting the path to the executables from the command line
idxdepth=$1

# Creating a variable for the reference genome
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Running idxdepth in parallel on 8 cores for all 102 samples
# DEPENDENCY : utilities/all_lines.txt
# DEPENDENCY : Illumina reads aligned with BWA
cat ../utilities/all_lines.txt | parallel -j 8 "$idxdepth -b ../illumina_data/aligned_reads/{}/{}_all.sort.bam -r $refgenome > idxdepth/{}.idxdepth"

# We take all the idxdepth files in the directory and compute the mean depth over all
# chromosomes for every sample. Has been called as "./gather_idxdepth.sh > idxdepth.txt"
for i in $(cat ../utilities/all_lines.txt)
do
	depth=$(grep "\"depth\"" idxdepth/${i}.idxdepth | egrep -o "[0-9]+\.[0-9]+" | awk '{x+=$1} END {print x / NR}')
	printf "$i $depth\n" >> idxdepth/idxdepth.txt
done


# From these data, we prepare one manifest file per sample in the manifest_files directory
for i in $(cat ../utilities/all_lines.txt)
do
    # Getting all the values to output to the file
    # WARNING : not sure if the relative path is going to work fine here
    filepath=../../illumina_data/aligned_reads/${i}/${i}_all.sort.bam
    depth=$(grep $i idxdepth/idxdepth.txt | awk '{print $2}')
    # DEPENDENCY : utilities/read_lengths.txt
    length=$(grep $i ../utilities/read_lengths.txt | awk '{print $3}')
    
    # Outputting the values
    printf "id\tpath\tdepth\tread length\n$i\t$filepath\t$depth\t$length\n" > manifest_files/${i}_manifest.txt
done

