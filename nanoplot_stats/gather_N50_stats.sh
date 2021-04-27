#!/bin/bash
printf "sample N50\n"

for i in $(cat /home/malem420/analyse_nanopore/genotyped_lines.txt)
do
	# Printing the name of the sample
	printf "$i " 

	# Printing the read N50
	N50=$(egrep --no-filename "Read length N50" /home/malem420/analyse_nanopore/${i}/guppy_4.0.11/nanoplots/sorted_bam_Gmaxv4/NanoStats.txt | egrep -o "[0-9,\.]*$")
	N50=${N50//,/}
	printf "$N50\n"
done

