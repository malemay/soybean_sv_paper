#!/bin/bash

# Collecting the data on SV counts per SVTYPE and SVLEN in order to generate Table S3

# Printing the header of the file
printf "sample smalldel mediumdel largedel hugedel smallins mediumins largeins hugeins smalldup mediumdup largedup hugedup smallinv mediuminv largeinv hugeinv\n" > table_s3_data.txt

# DEPENDENCY : utilities/line_ids.txt
for i in $(tail -n+2 ../utilities/line_ids.txt | cut -f1)
do
	# DEPENDENCY : scripts/count_svtypes_svsizes.awk
	# DEPENDENCY : filtered, refined and normalized SVs called by Sniffles
	sample=../nanopore_sv_calling/${i}_normalized_ids.vcf
	printf "$i " >> table_s3_data.txt
	../scripts/count_svtypes_svsizes.awk $sample >> table_s3_data.txt
done

