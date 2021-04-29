#!/bin/bash

# This code extracts the metadata on the Soybase TE sequences so it can be used
# in analyses. Essentially, we want to format the lines of the fasta with the
# information about the sequences in a tabular format so it can be input
# and manipulated within R

# I want the output file to contain the following fields :
# name | class | subclass | order | superfamily | family | description

grep "^>" te_database/SoyBase_TE_Fasta.txt | 
	awk '
	BEGIN{printf "name\tclass\tsubclass\torder\tsuperfamily\tfamily\tdescription\n"}
	{name = gensub("^.*name=([[:graph:]]+).*$", "\\1", "g"); printf name "\t"}
	{class = gensub("^.*\\sClass=([[:graph:]]+).*$", "\\1", "g"); printf class "\t"}
	{subclass = gensub("^.*Sub_Class=([[:graph:]]+).*$", "\\1", "g"); printf subclass "\t"}
	{order = gensub("^.*Order=([[:graph:]]+).*$", "\\1", "g"); printf order "\t"}
	{superfamily = gensub("^.*Super_Family=([[:graph:]]+).*$", "\\1", "g"); printf superfamily "\t"}
	{family = gensub("^.*\\sFamily=([[:graph:]]+).*$", "\\1", "g"); printf family "\t"}
	{description = gensub("^.*Description=([[:graph:]]*).*$", "\\1", "g"); if(!description) {description = "NA"} printf description "\n"}' \
		> te_metadata.txt

