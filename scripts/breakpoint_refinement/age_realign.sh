#!/bin/bash

# This script executes the assembly of the Oxford Nanopore reads,
# followed by polishing and re-alignment to the reference genome
# using AGE. It is meant to be called by the R function call_age()
# which will build the right command for it to be called
# getopts is used to parse the options sent from R

# Parsing the command options with getopts
while getopts s:m:a:w:c:t:h:e:l:r:b:n:f:p:z:x:y:u:d: opt; do

	case $opt in
	# Variables for the executables
	s) samtools=${OPTARG}
		;;
	m) minimap2=${OPTARG}
		;;
	a) age=${OPTARG}
		;;
	w) wtdbg2=${OPTARG}
		;;
	c) wtpoa_cns=${OPTARG}
		;;
	
	# Properties of the SV
	t) svtype=${OPTARG}
		;;
	h) chr=${OPTARG}
		;;
	e) start=${OPTARG}
		;;
	l) svlen=${OPTARG}
		;;
	
	# Input files
	r) refgenome=${OPTARG}
		;;
	b) nanopore_bam=${OPTARG}
		;;

	# Output files
	n) reads_fasta=${OPTARG}
		;;
	f) contig_fasta=${OPTARG}
		;;
	p) polished_fasta=${OPTARG}
		;;
	z) reference_fasta=${OPTARG}
		;;
	x) age_file=${OPTARG}
		;;
	y) prefix=${OPTARG}
		;;
	# Ranges for read extraction
	u) reads_window=${OPTARG}
		;;
	d) reference_window=${OPTARG}
		;;
	esac
done

# Some parameters are hard-coded
genome_size=500000
assembly_S=1
assembly_e=2
assembly_L=1000
gap_opening=-1

# Setting end_pos depending on the SV type
if [[ $svtype = DEL ]]
then
	let end_pos=$start+$svlen
elif [[ $svtype = INS ]]
then
	let end_pos=$start
fi

# Querying the reads overlapping that region
let query_start=$start-$reads_window
let query_end=$end_pos+$reads_window
$samtools view -b $nanopore_bam ${chr}:${query_start}-${query_end} | $samtools fasta - > $reads_fasta

# Assembling those reads using wtdbg2 and deriving consensus
$wtdbg2 -t 1 -x ont -S $assembly_S -e $assembly_e -L $assembly_L -g $genome_size -i $reads_fasta -f -o $prefix
$wtpoa_cns -t 1 -i ${prefix}.ctg.lay.gz -fo ${contig_fasta}

# All intermediate files can now be removed
rm ${prefix}*

# Calculating the number of contigs and exiting if not > 0
# The R script should handle the situation gracefully
# by detecting from the age file that no alignment was done
n_contigs=$(grep -c '^>' $contig_fasta)
if [[ $n_contigs -lt 1 ]]
then
	exit 0
fi

# Mapping the ONT reads with minimap and using them to polish
$minimap2 -t1 -a -x map-ont -r2k $contig_fasta $reads_fasta | $samtools sort | $samtools view -F0x900 - | $wtpoa_cns -t1 -d $contig_fasta -i - -fo $polished_fasta

# Extracting the sequences downstream and upstream of the region of interest for alignment
let subject_start=$start-$reference_window
let subject_end=$end_pos+$reference_window
$samtools faidx $refgenome ${chr}:${subject_start}-${subject_end} > $reference_fasta

# Running LongAGE on those two sequences
$age -allpos -go=${gap_opening} -indel -both $reference_fasta $polished_fasta > $age_file

