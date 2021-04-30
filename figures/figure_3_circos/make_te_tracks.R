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
# DEPENDENCY : refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3
repeats <- import("../../refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3")

# We extract the repeats of interest from the class column
ref_gypsy <- repeats[grepl("Gypsy", repeats$class)]
ref_copia <- repeats[grepl("Copia", repeats$class)]
ref_ltr   <- c(ref_gypsy, ref_copia)
ref_dna   <- repeats[grepl("DNA", repeats$class)]

# Let us look at the number of each of these
length(ref_gypsy)
length(ref_copia)
length(ref_ltr)
length(ref_dna)


# Now we also read the polymorphic elements from the analysis we made earlier
# DEPENDENCY: te_analysis/polymorphic_tes.tsv
polymorphic <- read.table("../../te_analysis/polymorphic_tes.tsv", 
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
length(poly_copia)
length(poly_ltr)
length(poly_dna)


# Now we can read in the genome bins in order to compute the number of TEs that overlap them
# DEPENDENCY : figures/figure_3_circos/gmax4_3Mb_bins.RData
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

# UTPUT: figures/figure_3_circos/ref_copia.txt
#output_track(ref_copia, gmax4_3Mb_bins, "ref_copia.txt")
# UTPUT: figures/figure_3_circos/ref_gypsy.txt
#output_track(ref_gypsy, gmax4_3Mb_bins, "ref_gypsy.txt")
# OUTPUT: figures/figure_3_circos/ref_ltr.txt
output_track(ref_ltr, gmax4_3Mb_bins, "ref_ltr.txt")
# OUTPUT: figures/figure_3_circos/ref_dna.txt
output_track(ref_dna, gmax4_3Mb_bins, "ref_dna.txt")
# UTPUT: figures/figure_3_circos/poly_copia.txt
#output_track(poly_copia, gmax4_3Mb_bins, "poly_copia.txt")
# UTPUT: figures/figure_3_circos/poly_gypsy.txt
#output_track(poly_gypsy, gmax4_3Mb_bins, "poly_gypsy.txt")
# OUTPUT: figures/figure_3_circos/poly_ltr.txt
output_track(poly_ltr, gmax4_3Mb_bins, "poly_ltr.txt")
# OUTPUT: figures/figure_3_circos/poly_dna.txt
output_track(poly_dna, gmax4_3Mb_bins, "poly_dna.txt")

