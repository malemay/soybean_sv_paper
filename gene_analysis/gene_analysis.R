#!/usr/bin/Rscript

# Loading the required packages
library(ggplot2)
library(GenomicFeatures)
library(VariantAnnotation)

# We will test the overlap between the SVs and genic features

# First I read the filtered population-scale Paragraph dataset
# We only need the AF and SVTYPE
vcf <- readVcf("/home/malem420/analyse_nanopore/transposons/paragraph_filtered.vcf",
	       param = ScanVcfParam(info = c("AF", "SVTYPE"), geno = NA))
vcf_ranges <- rowRanges(vcf)
vcf_ranges$af <- as.numeric(info(vcf)$AF)
vcf_ranges$svtype <- info(vcf)$SVTYPE
vcf_ranges$REF <- NULL
vcf_ranges$ALT <- NULL
rm(vcf)

# Now we read the gene annotation
annotations <- makeTxDbFromGFF("/home/malem420/refgenome/Gmax_v4/Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3")
# From the TxDb object we extract the objects for which we will want to test overlap
genes <- genes(annotations)
cds <- cds(annotations)
upstream5kb <- promoters(genes, upstream = 5000, downstream = 0)

# Saving these objects as .RData objects for retrieval later on
save(genes, file = "genes.RData")
save(cds, file = "cds.RData")
save(upstream5kb, file = "upstream5kb.RData")

# We initialize an "overlap" column that indicates what features the SV overlaps
vcf_ranges$overlap <- "intergenic"

# Then we become increasingly specific as to the type of feature it overlaps
vcf_ranges$overlap[overlapsAny(vcf_ranges, upstream5kb)] <- "upstream5kb"
vcf_ranges$overlap[overlapsAny(vcf_ranges, genes)] <- "gene"
vcf_ranges$overlap[overlapsAny(vcf_ranges, cds)] <- "cds"

# Tabulating the results
table(vcf_ranges$overlap)
# 
#         cds        gene  intergenic upstream5kb 
#        2071        8318       29686       16820 

table(vcf_ranges$overlap, vcf_ranges$svtype)
#              
#                 DEL   INS
#   cds          1685   386
#   gene         4329  3989
#   intergenic  15793 13893
#   upstream5kb  8853  7967

# Plotting the results to see whether there are any tendencies in allele frequency depending on the features that are overlapped
overlap_data <- as.data.frame(vcf_ranges)

# Saving the overlap_data as an R object for retrieval later
write.table(overlap_data, file = "overlap_data.txt", col.names = TRUE, row.names = FALSE,
	    quote = FALSE, sep = "\t")


# We also want to generate a table of the number of bases represented by each of these regions
# For this, we need to get a representation of the reference genome
# Reading the fasta index of the reference genome and pruning the unanchored scaffolds
refgenome <- scanFaIndex("/home/malem420/refgenome/Gmax_v4/Gmax_508_v4.0_mit_chlp.fasta")
chrs <- paste0("Gm", sprintf("%0.2d", 1:20))
seqlevels(refgenome, pruning.mod = "coarse") <- chrs

# We initialize a numeric vector of length 4 that will hold the number of bases corresponding to each feature type
feature_widths <- numeric(4)
names(feature_widths) <- c("cds", "gene", "intergenic", "upstream5kb")

# "cds" is simply the sum of the widths of the cds object, which we restrict to the reference chromosomes
seqlevels(cds, pruning.mode = "coarse") <- chrs
seqinfo(cds) <- seqinfo(refgenome)
feature_widths["cds"] <- sum(width(reduce(cds, ignore.strand = TRUE)))

# "gene" is the genic sequences that are NOT cds
seqlevels(genes, pruning.mode = "coarse") <- chrs
seqinfo(genes) <- seqinfo(refgenome)
feature_widths["gene"] <- sum(width(setdiff(genes, cds, ignore.strand = TRUE)))

# "upstream5kb" are the upstream sequences that do not overlap with genes
seqlevels(upstream5kb, pruning.mode = "coarse") <- chrs
seqinfo(upstream5kb) <- seqinfo(refgenome)
upstream5kb <- trim(upstream5kb)
feature_widths["upstream5kb"] <- sum(width(setdiff(upstream5kb, genes, ignore.strand = TRUE)))

# "intergenic" are the sequences that overlap neither upstream5kb nor genes
feature_widths["intergenic"] <- sum(width(setdiff(refgenome, union(upstream5kb, genes, ignore.strand = TRUE), ignore.strand = TRUE)))

# Just a sanity check
sum(width(refgenome)) == sum(feature_widths)

# Saving the object to file
save(feature_widths, file = "feature_widths.RData")

