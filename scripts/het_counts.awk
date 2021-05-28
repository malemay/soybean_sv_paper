#!/usr/bin/gawk -f
!/^#/ {# First determine the number of alternative alleles
      x = split($5, alleles, ",")
     
      # This script does not support cases with more than 10 alternative alleles
      if (x >= 10) {
          printf $1 " " $2 " " $4 " " $5 " " ".:.:." "\n"
          next
      }
 
      # Initialize an array with all possible values
      for(i = 0; i <= x + 1; i++) {
          counts[i]=0
      }

      # Iterate over all genotypes
      for(i = 10; i <= NF; i++) {
          a1 = substr($i, 1, 1)
          a2 = substr($i, 3, 1)
          # Incrementing the corresponding value (last field is heterozygous genotypes)
          if(a1 == a2 && a1 != ".") {counts[a1]++} else if(a1 != a2) {counts[x + 1]++}
      }
      
      # Printing the marker information to file (chr, position, alleles)
      printf $1 " " $2 " " $4 " " $5 " "
      # Then printing each of the genotype counts 
      for(i = 0; i <= x; i++) {
       printf counts[i] ":"
      }
      # Printing the last value without a following colon and with a new line
      printf counts[x + 1] "\n" 
     }

