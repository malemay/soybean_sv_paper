#!/bin/bash

# Getting the path to the SOAP executable from the command line
soap=$1

# DEPENDENCY : utilities/all_lines.txt
lines=$(cat ../../utilities/all_lines.txt)

# Moving to the soap directory
cd soap_assembly

for i in $lines
do
	mkdir -p $i
	./make_config.sh $i
	cd $i
	$soap sparse_pregraph -p 1 -s ${i}_config.txt -K 49 -R -z 1100000000 -o $i 1>pregraph.log 2>pregraph.err
	$soap contig -p 1 -g $i -R 1>contig.log 2>contig.err
	cd ..
done

