#!/bin/bash

# Getting the paths to the executables from the command line
blastn=$1
makeblastdb=$2

# Extracting the sequences of the deletions and insertions from the query vcf
# to blast them on the transposon database

# DEPENDENCY : te_analysis/query_all.vcf
grep -v "^#" query_all.vcf | awk '
$4 !~ /N/ && $5 !~ /N/ {
	if((length($4) - length($5)) >= 100) {
		print ">" $3
		print $4
	} else if ((length($5) - length($4)) >= 100) {
		print ">" $3
		print $5
        }
}
' > query_svs.txt

# Preparing the BLAST database
# DEPENDENCY : te_analysis/te_database/SoyBase_TE_Fasta.txt
$makeblastdb -dbtype nucl -in te_database/SoyBase_TE_Fasta.txt

# Blasting the query sequences against the TE database
$blastn -query query_svs.txt -db te_database/SoyBase_TE_Fasta.txt \
	-outfmt "6 qseqid qlen sseqid slen length qstart qend sstart send bitscore evalue pident" \
	-num_alignments 10 > blast_svs.txt

