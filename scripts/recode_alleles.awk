#!/usr/bin/gawk -f

# The explicitly coded structural variants result in plink files that are way too big
# Therefore, here I will change all REF alleles to A and ALT alleles to T, just for the
# purpose of recoding for population genetics analyses

# Changing the output field separator to \t
BEGIN {OFS = "\t"}

# Header lines are printed as is
/#/ {print $0}

# for other lines, REF and ALT alleles are changed and then the line is printed
!/^#/ {
	$4="A"
	$5="T"
	print
}

