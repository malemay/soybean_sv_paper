#!/bin/bash

# Getting the paths to the executables from the command line
bgzip=$1
tabix=$2
bcftools=$3

# Creating a variable for the reference fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Creating a variable for the regions to keep from the indel files
regions=Gm01,Gm02,Gm03,Gm04,Gm05,Gm06,Gm07,Gm08,Gm09,Gm10,Gm11,Gm12,Gm13,Gm14,Gm15,Gm16,Gm17,Gm18,Gm19,Gm20

# Looping over all the samples
for i in $(cat ../../utilities/all_lines.txt)
do
	# Assigning the files to process to variables
	sv_file=${i}.converted.vcf
	indel_file=${i}.svaba.indel.vcf

	# Removing the contig headers from each of the files and changing the Number for PL in the header from . to G
	grep -v "^##contig" $sv_file    | sed 's/##FORMAT=<ID=PL,Number=\./##FORMAT=<ID=PL,Number=G/' > sv_tmp.vcf

	# For the indels, we also add the annotation for the SVTYPE
	# DEPENDENCY : illumina_sv_calling/svaba/annotate_svtype.awk
	grep -v "^##contig" $indel_file | sed 's/##FORMAT=<ID=PL,Number=\./##FORMAT=<ID=PL,Number=G/' | ./annotate_svtype.awk > indel_tmp.vcf

	# Adding the headers
	$bcftools reheader --fai ${refgenome}.fai sv_tmp.vcf > sv_header.vcf
	$bcftools reheader --fai ${refgenome}.fai indel_tmp.vcf > indel_header.vcf

	# Compressing indel_header.vcf so I can filter out the scaffolds and organellar genomes
	$bgzip indel_header.vcf
	$tabix indel_header.vcf.gz

	# Processing the result through bcftools norm, and keeping only genotypes and the SV annotation
	$bcftools norm -f $refgenome -Ou sv_header.vcf | $bcftools view -G -Ou - | $bcftools annotate -x "^INFO/SVTYPE" -Ov - > ${i}_sv_norm.vcf
	$bcftools view -G --regions $regions -Ou indel_header.vcf.gz | $bcftools norm -f $refgenome -Ou - | $bcftools annotate -x "^INFO/SVTYPE" -Ov - > ${i}_indel_norm.vcf

	# Removing the temporary files
	rm sv_tmp.vcf indel_tmp.vcf sv_header.vcf indel_header.vcf.gz indel_header.vcf.gz.tbi
done

# Zipping and indexing all the files with bgzip and tabix so we can use bcftools merge on them
ls *indel_norm.vcf | parallel -j1 "$bgzip {}"
ls *sv_norm.vcf | parallel -j1 "$bgzip {}"

ls *indel_norm.vcf.gz | parallel -j1 "$tabix {}"
ls *sv_norm.vcf.gz | parallel -j1 "$tabix {}"

# Merging all SvABA normalized files using bcftools merge, and also adjusting the annotations
# DEPENDENCY : illumina_sv_calling/svaba/ACO_header_line.txt
$bcftools merge -m none -Ou $(ls *sv_norm.vcf.gz) $(ls *indel_norm.vcf.gz) | $bcftools norm -d none -Ov - | \
	gawk 'BEGIN {OFS = "\t"} /^#/ {print $0} !/^#/ {print $0 ";ACO=svaba"}' | \
	$bcftools annotate --header-lines ACO_header_line.txt -Ov - > svaba_merged.vcf

# This code extracts the SVs >= 50 nucleotides from the svaba variants
# It also sets the ID to the name of the caller + SV type + line number
# DEPENDENCY : scripts/extract_svs_50.awk
../../scripts/extract_svs_50.awk svaba_merged.vcf | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$3 = NR ; print}' | \
	$bcftools annotate --set-id "%INFO/ACO\_%INFO/SVTYPE\_%ID" -Ov - > svaba_svs.vcf

