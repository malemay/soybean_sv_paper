#!/bin/bash

# Getting the path to the executables from the command line
bgzip=$1
tabix=$2
vg=$3

# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
ref=../../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Creating a symbolic link to the vcf file in the current directory, then compressing and indexing it
# DEPENDENCY : breakpoint_refinement_analysis/refined_svs/svmerged_clustered_sorted.vcf
ln -s ../svmerged_clustered_sorted.vcf
$bgzip svmerged_clustered_sorted.vcf
$tabix svmerged_clustered_sorted.vcf.gz

# Running the command with option -a to ensure the paths will be kept for genotyping
seq -w 1 20 | xargs -I {} echo Gm{} | parallel -j 2 "$vg construct -a -t1 -m32 -C -R {} -r $ref -v svmerged_clustered_sorted.vcf.gz --flat-alts > {}.vg"

# Then converting the graphs to PackedGraph format in order to use less memory
ls Gm??.vg | parallel -j 2 "$vg convert -p {} > {.}.vgp"

# Making sure the ids are consistent across the graphs
$vg ids -j $(ls Gm??.vgp)

# Indexing the graph
$vg index -L -p -k 16 -x graph.xg $(ls Gm??.vgp) 

# Pruning the graph
ls -1 Gm??.vgp | parallel -j 4 "$vg prune -p -t1 -r {} > {.}.pruned.vgp"

# Indexing the graph with the GCSA index
mkdir -p tmpdir
$vg index -p -t 1 -b tmpdir -g graph.gcsa $(ls Gm??.pruned.vgp)

# Looping over all the samples
# DEPENDENCY : breakpoint_refinement_analysis/illumina_ids.txt
for i in $(cat ../../illumina_ids.txt)
do

# Creating variables for the location of the input files
# DEPENDENCY : trimmed Illumina files
fq1=../../../illumina_data/trimmed_fastq/${i}/${i}_R1_trimmed.fastq.gz
fq2=../../../illumina_data/trimmed_fastq/${i}/${i}_R2_trimmed.fastq.gz
fq3=../../../illumina_data/trimmed_fastq/${i}/${i}_sing_trimmed.fastq.gz

# Mapping the paired reads
$vg map -t 8 -d graph -f $fq1 -f $fq2 > ${i}_paired.gam

# Mapping the reads that were left unpaired by trimming with bbduk
$vg map -t 8 -d graph -f $fq3 > ${i}_single.gam

# Packing the alignments
$vg pack -t 4 -Q 5 -x graph.xg -g ${i}_paired.gam -o ${i}_paired.pack
$vg pack -i ${i}_paired.pack -t 4 -Q 5 -x graph.xg -g ${i}_single.gam -o ${i}_all.pack

# Calling variants based on the variant-aware graph
$vg call -t 4 -v svmerged_clustered_sorted.vcf.gz -k ${i}_all.pack graph.xg > ${i}_calls.vcf

# Compressing the vcf with bgzip
$bgzip ${i}_calls.vcf
$tabix ${i}_calls.vcf.gz

done

