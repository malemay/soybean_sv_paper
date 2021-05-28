#!/prg/R/4.0/bin/R

# Loading the required packages
library(GenomicFeatures)
library(VariantAnnotation)
library(rtracklayer)

# First I read the filtered population-scale Paragraph dataset
# We only need to read the allele frequencies (AF) and SV type (SVTYPE)
# DEPENDENCY: sv_genotyping/combined_svs/combined_paragraph_filtered.vcf
vcf <- readVcf("../sv_genotyping/combined_svs/combined_paragraph_filtered.vcf",
	       param = ScanVcfParam(fixed = NA,
				    info = c("AF", "SVTYPE"),
				    geno = NA))
vcf_ranges <- rowRanges(vcf)
vcf_ranges$af <- as.numeric(info(vcf)$AF)
vcf_ranges$svtype <- info(vcf)$SVTYPE
vcf_ranges$REF <- NULL
rm(vcf)

# Now let us create a GRanges representation of the reference genome
# This representation only includes the 20 reference chromosomes
# and excludes the unanchored scaffolds

# Reading the fasta index of the reference genome and pruning the unanchored scaffolds
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome <- scanFaIndex("../refgenome/Gmax_508_v4.0_mit_chlp.fasta")
chrs <- paste0("Gm", sprintf("%0.2d", 1:20))
seqlevels(refgenome, pruning.mod = "coarse") <- chrs

# Now we either read the annotation from disk or generate it from the GFF3 file
if(!file.exists("genes.RData") || !file.exists("cds.RData") || !file.exists("upstream5kb.RData")) {
	# Generating the database from the GFF3 file
	# DEPENDENCY : refgenome/Gmax_508_Wm82.a4.v1.gene_exons.gff3
	annotations <- makeTxDbFromGFF("../refgenome/Gmax_508_Wm82.a4.v1.gene_exons.gff3")

	# From the TxDb object we extract the objects for which we will want to test overlap
	# We also restrict them to those found on the reference chromosomes
	genes <- genes(annotations)
	seqlevels(genes, pruning.mode = "coarse") <- chrs
	seqinfo(genes) <- seqinfo(refgenome)

	cds <- cds(annotations)
	seqlevels(cds, pruning.mode = "coarse") <- chrs
	seqinfo(cds) <- seqinfo(refgenome)

	upstream5kb <- promoters(genes, upstream = 5000, downstream = 0)
	seqlevels(upstream5kb, pruning.mode = "coarse") <- chrs
	seqinfo(upstream5kb) <- seqinfo(refgenome)
	# The trimming here makes sure that we only keep sequences that are within the bounds of the reference genome
	upstream5kb <- trim(upstream5kb)

	# Saving these objects as .RData objects for retrieval later on
	save(genes, file = "genes.RData")
	save(cds, file = "cds.RData")
	save(upstream5kb, file = "upstream5kb.RData")
} else {
	load("genes.RData")
	load("cds.RData")
	load("upstream5kb.RData")
}

# We initialize an "overlap" column that indicates what features the SV overlaps
vcf_ranges$overlap <- "intergenic"

# Then we become increasingly specific as to the type of feature it overlaps
vcf_ranges$overlap[overlapsAny(vcf_ranges, upstream5kb)] <- "upstream5kb"
vcf_ranges$overlap[overlapsAny(vcf_ranges, genes)] <- "gene"
vcf_ranges$overlap[overlapsAny(vcf_ranges, cds)] <- "cds"

# Coercing the overlap data to a data.frame for plotting with ggplot2
overlap_data <- as.data.frame(vcf_ranges)

# Saving the overlap_data as text file for retrieval from disk later on
# OUTPUT : gene_analysis/overlap_data.txt
write.table(overlap_data, file = "overlap_data.txt", 
	    col.names = TRUE, row.names = FALSE,
	    quote = FALSE, sep = "\t")

# We initialize a numeric vector of length 4 that will hold the number of bases corresponding to each feature type
feature_widths <- numeric(4)
names(feature_widths) <- c("cds", "gene", "intergenic", "upstream5kb")

# "cds" is simply the sum of the widths of the cds object, which we restrict to the reference chromosomes
cds_ranges <- reduce(cds, ignore.strand = TRUE)
feature_widths["cds"] <- sum(width(cds_ranges))

# "gene" is the genic sequences that are NOT cds
gene_ranges <- setdiff(genes, cds, ignore.strand = TRUE)
feature_widths["gene"] <- sum(width(gene_ranges))

# "upstream5kb" are the upstream sequences that do not overlap with genes
upstream5kb_ranges <- setdiff(upstream5kb, genes, ignore.strand = TRUE)
feature_widths["upstream5kb"] <- sum(width(upstream5kb_ranges))

# "intergenic" are the sequences that overlap neither upstream5kb nor genes
intergenic_ranges <- setdiff(refgenome, union(upstream5kb, genes, ignore.strand = TRUE), ignore.strand = TRUE)
feature_widths["intergenic"] <- sum(width(intergenic_ranges))

# Just a sanity check that the sums of the computed widths corresponds to the size of the genome
stopifnot(sum(width(refgenome)) == sum(feature_widths))

# Saving the object to file
# OUTPUT : gene_analysis/feature_widths.RData
save(feature_widths, file = "feature_widths.RData")

# Let us restrict our attention to the non-repeated regions
# We can simply read the bed file that we generated before for the benchmarks in non-repeat regions
# DEPENDENCY : refgenome/repeat_regions/non_repeated_regions.bed
non_repeats <- import("../refgenome/repeat_regions/non_repeated_regions.bed")

# We can find the intersection between each set of genic feature ranges and non_repeat ranges
# in order to find the regions of each that are in non-repeat regions
cds_norepeat <- intersect(cds_ranges, non_repeats, ignore.strand = TRUE)
gene_norepeat <- intersect(gene_ranges, non_repeats, ignore.strand = TRUE)
upstream5kb_norepeat <- intersect(upstream5kb_ranges, non_repeats, ignore.strand = TRUE)
intergenic_norepeat <- intersect(intergenic_ranges, non_repeats, ignore.strand = TRUE)

# Filling a named numeric vector with the total length of every genic feature in non-repeat regions
norepeat_widths <- numeric(4)
names(norepeat_widths) <- c("cds", "gene", "intergenic", "upstream5kb")

norepeat_widths["cds"] <- sum(width(cds_norepeat))
norepeat_widths["gene"] <- sum(width(gene_norepeat))
norepeat_widths["upstream5kb"] <- sum(width(upstream5kb_norepeat))
norepeat_widths["intergenic"] <- sum(width(intergenic_norepeat))

# Just a sanity check
stopifnot(sum(norepeat_widths) == sum(width(non_repeats)))

# Saving this to file for use in reporting results
# OUTPUT : gene_analysis/norepeat_widths.RData
save(norepeat_widths, file = "norepeat_widths.RData")

# Now we can restrict the SVs to those that overlap non-repeat regions
norepeat_vcf <- subsetByOverlaps(vcf_ranges, non_repeats)

# We re-initialize the "overlap" column that indicates what features the SV overlaps
norepeat_vcf$overlap <- "intergenic"

# Then we become increasingly specific as to the type of feature it overlaps
norepeat_vcf$overlap[overlapsAny(norepeat_vcf, upstream5kb_norepeat)] <- "upstream5kb"
norepeat_vcf$overlap[overlapsAny(norepeat_vcf, gene_norepeat)] <- "gene"
norepeat_vcf$overlap[overlapsAny(norepeat_vcf, cds_norepeat)] <- "cds"

# Coercing the non-repeat overlap data to a data.frame
overlap_norepeat <- as.data.frame(norepeat_vcf)

# Also saving this table to file for reporting the results
# OUTPUT : gene_analysis/overlap_norepeat.txt
write.table(overlap_norepeat, file = "overlap_norepeat.txt", 
	    col.names = TRUE, row.names = FALSE,
	    quote = FALSE, sep = "\t")

