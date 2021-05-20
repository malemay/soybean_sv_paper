#!/bin/bash

# Getting the path to the bcftools executable from the command line
bcftools=$1

# We sort the vcfs and remove all the annotations that we do not need anymore
for i in illumina nanopore
do
	$bcftools sort -m 1000.0M -Ou ${i}_merged.vcf | \
	       	$bcftools annotate -x "^INFO/ACO,INFO/InitialIDs,INFO/SVTYPE,INFO/ClusterIDs" -Ov - > ${i}_merged_sorted.vcf
done

