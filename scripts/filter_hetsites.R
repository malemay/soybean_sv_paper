#!/usr/bin/Rscript
hetrates <- read.table("hetrates.txt", sep = " ", colClasses = c("character", "numeric", "NULL", "NULL", "NULL"))
hetstats <- read.table("wgs_hetstats.txt", sep = " ", header = TRUE, colClasses = c("numeric", "logical", "numeric"))

if(nrow(hetrates) != nrow(hetstats)) stop("ERROR: number of rows not equal")

hetrates <- hetrates[hetstats$HETRATE <= 0.1 & hetstats$HOMPASS, ]

write.table(hetrates, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

