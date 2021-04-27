#!/usr/bin/gawk -f

# This script takes a vcf file as input and returns a line with
# the number of SVs in a given size class for each SV type,
# according to the following fields :
# smalldel mediumdel largedel hugedel smallins mediumins largeins hugeins smalldup mediumdup largedup hugedup smallinv mediuminv largeinv hugeinv
# Where :
# 	small denotes variants [50-100 bp[
#	medium denotes variants [100-1,000bp[
#	large denotes variants [1,000-10,000bp[
#	huge denotes variants >= 10,000bp

# Initializing all values to 0, as a safety measure
BEGIN {smalldel = 0; mediumdel = 0; largedel = 0; hugedel = 0;
	smallins = 0; mediumins = 0; largeins = 0; hugeins = 0;
	smalldup = 0; mediumdup = 0; largedup = 0; hugedup = 0;
	smallinv = 0; mediuminv = 0; largeinv = 0; hugeinv = 0}

# Computing a generic SV length from the lengths of the REf and ALT alleles
{svlen = length($5) - length($4)}

# If it is a deletion
/SVTYPE=DEL/ && svlen <= -50 && svlen > -100 {smalldel+=1}
/SVTYPE=DEL/ && svlen <= -100 && svlen > -1000 {mediumdel+=1}
/SVTYPE=DEL/ && svlen <= -1000 && svlen > -10000 {largedel+=1}
/SVTYPE=DEL/ && svlen <= -10000 {hugedel+=1}

# If it is an insertion
/SVTYPE=INS/ && svlen >= 50 && svlen < 100 {smallins+=1}
/SVTYPE=INS/ && svlen >= 100 && svlen < 1000 {mediumins+=1}
/SVTYPE=INS/ && svlen >= 1000 && svlen < 10000 {largeins+=1}
/SVTYPE=INS/ && svlen >= 10000 {hugeins+=1}

# If it is a duplication
/SVTYPE=DUP/ && svlen >= 50 && svlen < 100 {smalldup+=1}
/SVTYPE=DUP/ && svlen >= 100 && svlen < 1000 {mediumdup+=1}
/SVTYPE=DUP/ && svlen >= 1000 && svlen < 10000 {largedup+=1}
/SVTYPE=DUP/ && svlen >= 10000 {hugedup+=1}

# If it is an inversion, then we get the size from the reference length
/SVTYPE=INV/ && length($4) >= 50 && length($4) < 100 {smallinv+=1}
/SVTYPE=INV/ && length($4) >= 100 && length($4) < 1000 {mediuminv+=1}
/SVTYPE=INV/ && length($4) >= 1000 && length($4) < 10000 {largeinv+=1}
/SVTYPE=INV/ && length($4) >= 10000 {hugeinv+=1}


END{printf smalldel " " mediumdel " " largedel " " hugedel " " smallins " " mediumins " " largeins " " hugeins " ";
	printf smalldup " " mediumdup " " largedup " " hugedup " " smallinv " " mediuminv " " largeinv " " hugeinv "\n"}

