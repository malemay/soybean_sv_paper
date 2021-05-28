#!/bin/bash

# Assigning the bash argument to a variable
n_ind=$1
# gawk is invoked with the bash argument bearing name n_ind
gawk -v n_ind="$n_ind" '

### - BEGIN STATEMENT
# The field separator is : for the genotype counts ; adding column names
BEGIN {FS = ":" ; printf "HETRATE HOMPASS MISSING\n"}

### - LINE PROCESSING
{
 # Initializing a value for the row sum
 rowsum = 0
 # This will hold the number of homozygous types that > heterozygotes
 hom_flag = 0

 # If the value is missing for this record then only missing values are output
 if($1 == ".") {
     printf "NA NA NA" "\n"
     next
 }

 # Iterating over each genotype count
 for(i = 1; i <= NF; i++) {
     # Incrementing the row sum
     rowsum+=$i
     # Incrementing the hom_flag if higher than the number of heterozygotes
     if ($i > $NF) hom_flag++
 }

 ### --- OUTPUT
 # Print the heterozygosity rate
 if (rowsum > 0) printf $NF / rowsum " "; else printf "NA "
 # Print the homozygous filter flag
 if (hom_flag >= 2) printf "TRUE " ; else printf "FALSE "
 # Print the number of iussing genotypes
 printf n_ind - rowsum
 # Printing a new line before the next iteration
 printf "\n"
}
'

