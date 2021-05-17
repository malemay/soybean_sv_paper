#!/bin/bash

# Getting the path to the AsmVar executable from the command line
asmvar=$1

# Creating a variable to hold the names of the 20 chromosomes
chrlist=$(seq -w 1 20 | xargs -I {} echo -n "Gm{} ")

# Moving to the asmvar_calling directory
cd asmvar_calling

# DEPENDENCY : utilities/all_lines.txt
for j in $(cat ../../../utilities/all_lines.txt)
do
 
  mkdir -p $j
  cd $j

  for i in $chrlist
  do
	  # DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
	  # DEPENDENCY : LAST alignment files
	  # DEPENDENCY : SOAPdenovo2 assemblies
	  $asmvar -s $j -i ../../last_alignment/${j}_soapdenovo2_fm.maf \
		  -t ../../../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
		  -q ../../soap_assembly/${j}/${j}.contig -o asmvar_results_${i} -r $i > ${j}_${i}.age
  done

  cd ..
done

