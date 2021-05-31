#!/prg/R/4.0/bin/Rscript

# File initially created on April 16, 2021

# This code contains a permutation test to see whether the proportions
# of insertions and deletions overlapping various genic features in our
# dataset differs from what would be randomly expected

# The randomization will be performed as follows:
#
# SV start positions will be randomly shuffled within 100-kb
# bins while keeping their width constant. The reason for
# performing the randomization within 100-kb bins is that the
# genome has local peculiarities that may be difficult to explicitly
# account for. By reshuffling the SV positions within those bins, we
# preserve the local SV distribution while allowing the positions of
# the SVs to change relative to genes.

# Loading the required libraries
library(GenomicRanges)
library(Rsamtools)
library(VariantAnnotation)
library(rtracklayer)

# This function takes the GRanges object of SVs and splits it into
# a GRangesList object with one element for each bin of the
# genome
split_ranges <- function(sv_ranges, refgenome, binwidth) {
	# Generating a tiling of the genome from the reference genome
	refbins <- tileGenome(seqinfo(refgenome), tilewidth = binwidth, cut.last.tile.in.chrom = TRUE)

	# Getting the first overlap within refbins of each element of sv_ranges
	overlaps <- findOverlaps(sv_ranges, refbins, select = "first")
	if(any(is.na(overlaps))) stop("Error while finding matches in genome tiles")

	# Splitting the sv_ranges object according to their matches within the tiled genome
	output <- split(sv_ranges, overlaps)

	# Setting the names of the output object to be of the format "chr:start-end"
	indices <- as.numeric(names(output))
	names(output) <- paste0(seqnames(refbins[indices]), ":", start(refbins[indices]), "-", end(refbins[indices]))
	output
}

# This function shuffles the positions of the input GRangesList element-wise and returns a GRanges object
shuffle_positions <- function(split_ranges) {
	# Extracting the bins starting and ending positions
	bin_names  <- sub("Gm[0-9]{2}:", "", names(split_ranges))
	bin_starts <- as.numeric(sapply(strsplit(bin_names, "-"), `[`, 1))
	bin_starts <- rep(bin_starts, times = lengths(split_ranges))
	bin_ends   <- as.numeric(sapply(strsplit(bin_names, "-"), `[`, 2))
	bin_ends <- rep(bin_ends, times = lengths(split_ranges))

	# Unlisting the split_ranges as a single GRanges object
	split_ranges <- unlist(split_ranges)

	# Re-setting the starts and widths
	widths <- width(split_ranges)
	starts <- round(runif(length(split_ranges), bin_starts, bin_ends))
	ranges(split_ranges) <- IRanges(start = starts, width = widths)

	split_ranges
}

# A function that takes a GRanges of SVs as well as GRanges of cds, genes, and upstream regions,
# and returns two numeric vectors with the proportions of deletions and insertions, respectively,
# overlapping various genic features
check_overlap <- function(sv_ranges, cds_ranges, gene_ranges, upstream_ranges, norepeats = NULL) {

	# We subset the vcf ranges if there is a norepeats GRangse object specified
	# We keep only those SVs that have an overlap to the norepeats GRanges
	if(!is.null(norepeats)) sv_ranges <- subsetByOverlaps(sv_ranges, norepeats)

	# We initialize an "overlap" column that indicates what features the SV overlaps
	sv_ranges$overlap <- "intergenic"

	# Then we become increasingly specific as to the type of feature it overlaps
	sv_ranges$overlap[overlapsAny(sv_ranges, upstream_ranges)] <- "upstream5kb"
	sv_ranges$overlap[overlapsAny(sv_ranges, gene_ranges)] <- "gene"
	sv_ranges$overlap[overlapsAny(sv_ranges, cds_ranges)] <- "cds"

	# Computing the number of deletions and insertions overlapping various features
	counts <- table(sv_ranges$overlap, sv_ranges$svtype)
	del_props <- counts[, "DEL"] / sum(counts[, "DEL"])
	ins_props <- counts[, "INS"] / sum(counts[, "INS"])

	list(del = del_props, ins = ins_props)
}

# Creating a function that replicates the analysis a given number of time,
# given some bin width
sv_permute <- function(sv_ranges, cds_ranges, gene_ranges, upstream_ranges, niter, binsize, norepeats = NULL) {

	sv_ranges_split <- split_ranges(sv_ranges, refgenome, binsize)

	# Preparing matrices to hold the output
	del <- matrix(nrow = niter, ncol = 4)
	ins <- matrix(nrow = niter, ncol = 4)

	for(i in 1:niter) {
		i_svs <- shuffle_positions(sv_ranges_split)
		sv_props <- check_overlap(i_svs, cds_ranges, gene_ranges, upstream_ranges, norepeats)
		del[i, ] <- sv_props$del
		ins[i, ] <- sv_props$ins
	}

	# Adding column names
	colnames(del) <- names(sv_props$del)
	colnames(ins) <- names(sv_props$ins)

	list(del = del, ins = ins)
}

# Now that we have all the functions, we can perform the analysis

# Creating a GRanges representation of the reference genome with only chromsomes Gm01 to Gm20
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome <- scanFaIndex("../refgenome/Gmax_508_v4.0_mit_chlp.fasta")
chrs <- paste0("Gm", sprintf("%0.2d", 1:20))
seqlevels(refgenome, pruning.mod = "coarse") <- chrs

# Reading the SV GRanges from the vcf file
# DEPENDENCY : sv_genotyping/combined_svs/combined_paragraph_filtered.vcf
vcf <- readVcf("../sv_genotyping/combined_svs/combined_paragraph_filtered.vcf",
	       param = ScanVcfParam(fixed = NA,
				    info = c("AF", "SVTYPE"),
				    geno = NA))
vcf_ranges <- rowRanges(vcf)
vcf_ranges$af <- as.numeric(info(vcf)$AF)
vcf_ranges$svtype <- info(vcf)$SVTYPE
vcf_ranges$REF <- NULL
rm(vcf)

# Reading the genes, cds and upstream5kb GRanges from file 
# DEPENDENCY : gene_analysis/genes.RData
load("genes.RData")
# DEPENDENCY : gene_analysis/cds.RData
load("cds.RData")
# DEPENDENCY : gene_analysis/upstream5kb.RData
load("upstream5kb.RData")

# First testing the functions to make sure that everything works as intended
#### ----- Testing the split_ranges function
split_vcf <- split_ranges(sv_ranges = vcf_ranges, refgenome = refgenome, binwidth = 500000)
# Length of the input preserved
stopifnot(length(unlist(split_vcf)) == length(vcf_ranges))
# seqnames preserved
stopifnot(all(seqnames(unlist(split_vcf)) == seqnames(vcf_ranges)))
# start positions do fall within the defined ranges
stopifnot(all(unlist(start(split_vcf)) >= as.numeric(rep(sapply(strsplit(sub(".*:", "", names(split_vcf)), "-"), `[`, 1), lengths(split_vcf)))))
stopifnot(all(unlist(start(split_vcf)) <= as.numeric(rep(sapply(strsplit(sub(".*:", "", names(split_vcf)), "-"), `[`, 2), lengths(split_vcf)))))

#### ----- Testing the shuffle_positions function
shuffle_test <- shuffle_positions(split_vcf)
# Length of the input preserved
stopifnot(length(vcf_ranges) == length(shuffle_test))
# Width of the variants preserved
stopifnot(all(width(vcf_ranges) == width(shuffle_test)))
# SVTYPE preserved
stopifnot(all(vcf_ranges$svtype == shuffle_test$svtype))
# Distribution of the differences in starting position
if(interactive()) hist(start(vcf_ranges) - start(shuffle_test))
stopifnot(all(abs(start(vcf_ranges) - start(shuffle_test)) <= 500000))

#### ----- Testing the check_overlap function
# Comparing the output to the one already computed for the real data
if(interactive()) check_overlap(vcf_ranges, cds, genes, upstream5kb)

#### ----- Testing the sv_permute function at small scale
if(interactive()) sv_permute(vcf_ranges, cds, genes, upstream5kb, 10, 500000)
# Everything seems to be working fine

# Now launching the analysis 
permutation_all_100kb <- sv_permute(vcf_ranges, cds, genes, upstream5kb, 5000, 100000)

# Saving these objects to file
# OUTPUT : gene_analysis/permutation_all_100kb.RData
save(permutation_all_100kb, file = "permutation_all_100kb.RData")

