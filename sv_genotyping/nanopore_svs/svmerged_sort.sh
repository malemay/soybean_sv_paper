#!/bin/bash

# Getting the path to the bcftools executable from the command line
bcftools=$1

# Sort vcf entries
$bcftools sort -m 500.0M -Ov svmerged.clustered.vcf > svmerged_clustered_sorted.vcf

