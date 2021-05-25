#!/bin/bash
printf "sample N50\n"

# DEPENDENCY : utilities/line_ids.txt
for i in $(tail -n +2 ../utilities/line_ids.txt | cut -f 1)
do
	# Printing the name of the sample
	printf "$i " 

	# Printing the read N50
	# DEPENDENCY : NanoPlot stats on the sorted bam files
	N50=$(egrep --no-filename "Read length N50" ${i}/NanoStats.txt | egrep -o "[0-9,\.]*$")
	N50=${N50//,/}
	printf "$N50\n"
done

