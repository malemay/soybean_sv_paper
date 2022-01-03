#!/bin/bash

# Getting the path to the executables from the command line
esearch=$1
efetch=$2

# Getting the metadata for the Illumina data from the SRA using BioProject #PRJNA356132
$esearch -db sra -query "PRJNA356132[GPRJ]" | $efetch -format runinfo | sort | grep -v ^$ | uniq > sra_metadata.csv

