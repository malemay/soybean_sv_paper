#!/usr/bin/Rscript

# This code creates the data tracks for plotting the number of transposable
# elements per bin using Circos

# Eight tracks need to be prepared
# - Number of Gypsy elements per bin in the reference genome
# - Number of Copia elements per bin in the reference genome
# - Number of LTR retrotransposons (Gypsy + Copia) per bin in the reference genome
# - Number of DNA transposable elements per bin in the reference genome
# - Number of polymorphic Gypsy elements per bin 
# - Number of polymorphic Copia elements per bin 
# - Number of polymorphic LTR retrotransposons (Gypsy + Copia) per bin
# - Number of polymorphic DNA transposable elements per bin 

# Loading the required packages
library(rtracklayer)
library(GenomicRanges)

# First we will read the repeat annotation from the repeat GFF file
repeats <- import("/home/malem420/refgenome/Gmax_v4/Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3")

# We extract the repeats of interest from the class column
ref_gypsy <- repeats[grepl("Gypsy", repeats$class)]
ref_copia <- repeats[grepl("Copia", repeats$class)]
ref_ltr   <- c(ref_gypsy, ref_copia)
ref_dna   <- repeats[grepl("DNA", repeats$class)]

# Let us look at the number of each of these
length(ref_gypsy)
# [1] 74930
length(ref_copia)
# [1] 96745
length(ref_ltr)
# [1] 171675
length(ref_dna)
# [1] 99866


# Now we also read the polymorphic elements from the analysis we made earlier
polymorphic <- read.table("/home/malem420/analyse_nanopore/transposons/polymorphic_tes.tsv", 
			  header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# Then we coerce it as a GRanges object
polymorphic <- makeGRangesFromDataFrame(polymorphic, keep.extra.columns = TRUE)
# We now extract the groups of interest as individual objects
poly_gypsy <- polymorphic[polymorphic$te_group == "Gypsy"]
poly_copia <- polymorphic[polymorphic$te_group == "Copia"]
poly_ltr   <- c(poly_gypsy, poly_copia)
poly_dna   <- polymorphic[polymorphic$te_group == "DNA"]

# Now let us have a look at their numbers
length(poly_gypsy)
# [1] 1622
length(poly_copia)
# [1] 2375
length(poly_ltr)
# [1] 3997
length(poly_dna)
# [1] 397


# Now we can read in the genome bins in order to compute the number of TEs that overlap them
load("gmax4_3Mb_bins.RData")

# Let us write a function that takes a GRanges object as input, some bins to count over
# and the name of a file to which the Circos data track will be output
output_track <- function(x, bins, output_file) {
	bins$count <- countOverlaps(bins, x)
	bins <- as.data.frame(bins)
	bins <- bins[, c("seqnames", "start", "end", "count")]
	write.table(bins, file = output_file, col.names = FALSE, row.names = FALSE,
		    quote = FALSE, sep = "\t")
}

output_track(ref_copia, gmax4_3Mb_bins, "ref_copia.txt")
output_track(ref_gypsy, gmax4_3Mb_bins, "ref_gypsy.txt")
output_track(ref_ltr, gmax4_3Mb_bins, "ref_ltr.txt")
output_track(ref_dna, gmax4_3Mb_bins, "ref_dna.txt")
output_track(poly_copia, gmax4_3Mb_bins, "poly_copia.txt")
output_track(poly_gypsy, gmax4_3Mb_bins, "poly_gypsy.txt")
output_track(poly_ltr, gmax4_3Mb_bins, "poly_ltr.txt")
output_track(poly_dna, gmax4_3Mb_bins, "poly_dna.txt")

