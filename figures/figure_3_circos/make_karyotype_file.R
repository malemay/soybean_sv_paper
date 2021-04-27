#!/usr/bin/Rscript

# This code prepares a karyotype file that specifies the span
# of the soybean reference chromosomes Gm01 to Gm20 and attributes
# default colors to them for plotting with Circos

# Loading the Rsamtools library for reading in the fasta index
library(Rsamtools)

# Reading in the index
fai <- scanFaIndex("/home/malem420/refgenome/Gmax_v4/Gmax_508_v4.0_mit_chlp.fasta")

# Removing organelles and unanchored genomes
fai <- fai[grepl("^Gm[0-9]{2}$", seqnames(fai))]

# Formatting the data for output to the Gmax_karyotype.txt file
fai <- as.data.frame(fai)
fai$seqnames <- as.character(fai$seqnames)
fai <- fai[, c("seqnames", "start", "end")]

fai$start <- fai$start - 1

names(fai) <- c("ID", "START", "END")

fai$chr <- "chr"
fai$sep <- "-"
fai$LABEL <- fai$ID
fai$COLOR <- paste0("chr", 1:20)

fai <- fai[, c("chr", "sep", "ID", "LABEL", "START", "END", "COLOR")]

# Writing to file
write.table(fai, file = "Gmax_karyotype.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

