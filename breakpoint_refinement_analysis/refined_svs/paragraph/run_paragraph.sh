#!/bin/bash

# Getting the path to the executables from the command line
bcftools=$1
bgzip=$2
tabix=$3
multigrmpy=$4

# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Creating a symbolic link to the input vcf file
# DEPENDENCY : breakpoint_refinement_analysis/refined_svs/svmerged_clustered_sorted.vcf
ln -s ../svmerged_clustered_sorted.vcf all_svs.vcf

# Paragraph won't accept variants that are less than a read length away
# from chromosome ends. This line will therefore extract the positions 1000 and END - 1000 for each chromsome
# in a tab-separated format. I chose 1000 bp instead of 100 because some indel sites were kept even though they started before 100
gawk '/^##contig=<ID=Gm..,/ {match($0, /Gm../, chr) ;  match($0, /[0-9]{6,}/, N) ; print chr[0] "\t" 1000 "\t" N[0] - 1000}' all_svs > regions_file.txt

# Compressing and indexing the file, otherwise the next command won't work
$bgzip all_svs.vcf
$tabix all_svs.vcf.gz
$bcftools view --regions-file regions_file.txt -Ov all_svs.vcf.gz > all_svs_truncated.vcf

# Padding the variants as required by Paragraph
# There might be an issue with the path to the reference genome in addMissingPaddingGmax4.py . Adjust if needed
# DEPENDENCY : scripts/addMissingPaddingGmax4.py
python2 ../../../scripts/addMissingPaddingGmax4.py all_svs_truncated.vcf > all_svs_padded.vcf 

# Removing the annotation columns that cause bugs or warnings when running paragraph
$bcftools annotate -x INFO/REF_strand,INFO/END,INFO/SUPTYPE -Ov all_svs_padded.vcf > all_svs_padded_ok.vcf

# For each iteration, we extract the -M option from the manifest file and then launch multigrmpy
# DEPENDENCY : breakpoint_refinement_analysis/illumina_ids.txt
# DEPENDENCY : manifest files
# DEPENDENCY : Illumina alignment files
cat ../../illumina_ids.txt | \
	parallel -j16 ' m_opt=$(grep {} ../../../sv_genotyping/manifest_files/{}_manifest.txt | awk '\''{print int($3 * 20)}'\'') ; \
	python3 $multigrmpy -t1 -M $m_opt -i all_svs_padded_ok.vcf -m ../../../sv_genotyping/manifest_files/{}_manifest.txt -r $refgenome --scratch-dir tmpdir/ -o {}_results ;'

