#!/bin/bash

# Collecting the data on SV counts per SVTYPE and SVLEN in order to generate Table S3

# Printing the header of the file
printf "sample smalldel mediumdel largedel hugedel smallins mediumins largeins hugeins smalldup mediumdup largedup hugedup smallinv mediuminv largeinv hugeinv\n" > table_s3_data.txt

# DEPENDENCY : genotyped_lines.txt
for i in $(cat /home/malem420/analyse_nanopore/genotyped_lines.txt)
do
	# DEPENDENCY : ../scripts/count_svtypes_svsizes.awk
	# DEPENDENCY : ${i}_sniffles_normalized_ids.vcf
	sample=/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_nanopore/genotyping/${i}_normalized_ids.vcf
	printf "$i " >> table_s3_data.txt
	../scripts/count_svtypes_svsizes.awk $sample >> table_s3_data.txt
done

