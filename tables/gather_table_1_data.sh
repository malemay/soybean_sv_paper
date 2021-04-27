#!/bin/bash

# Extracting metadata from the vcf files of each of the SV calling programs
#  in order to make table 1, which will summarize the number of variants of
#  each SV type and size class for each program and for the merged data set
#  This table is valid for the programs calling SVs from Illumina data.

# Looping over the files that contain all the SVs from every tool
# DEPENDENCY : asmvar_svs.vcf
# DEPENDENCY : manta_svs.vcf
# DEPENDENCY : smoove_svs.vcf
# DEPENDENCY : svaba_svs.vcf
# DEPENDENCY : svmerged.clustered.vcf
asmvar=/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_genotyping/asmvar_variants/asmvar_svs.vcf
manta=/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_genotyping/manta_variants/manta_svs.vcf
smoove=/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_genotyping/smoove_variants/smoove_svs.vcf
svaba=/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_genotyping/svaba_variants/svaba_svs.vcf
all=/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_genotyping/svmerged_variants/svmerged.clustered.vcf

# Printing the header of the file
printf "program smalldel mediumdel largedel hugedel smallins mediumins largeins hugeins smalldup mediumdup largedup hugedup smallinv mediuminv largeinv hugeinv\n" > table_1_data.txt

# DEPENDENCY : genotyped_lines.txt
for i in $asmvar $manta $smoove $svaba $all
do
	# DEPENDENCY : ../scripts/count_svtypes_svsizes.awk
	printf "$i " >> table_1_data.txt
	../scripts/count_svtypes_svsizes.awk $i >> table_1_data.txt
done

# Then reprocessing the file to change the values of the first column
gawk 'BEGIN {OFS = ","} 
	/program/ {$1 = $1}
	/asmvar/ {$1 = "asmvar"} 
	/manta/ {$1 = "manta"} 
	/smoove/ {$1 = "smoove"} 
	/svaba/ {$1 = "svaba"} 
	/svmerged/ {$1 = "merged"} 
	{print}' table_1_data.txt > table_1.csv

