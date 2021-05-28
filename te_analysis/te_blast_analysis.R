#!/prg/R/4.0/bin/Rscript

# Loading the required packages
library(VariantAnnotation)
library(rtracklayer)

# Reading the VCF containing the information on the queried SVs
# DEPENDENCY : te_analysis/query_all.vcf
vcf <- readVcf("query_all.vcf", param = ScanVcfParam(info = c("SVTYPE", "AF"), geno = NA))

# Reading the results of the blastn query
# DEPENDENCY : te_analysis/blast_svs.txt
blastn_results <- read.table("blast_svs.txt", header = FALSE, stringsAsFactors = FALSE,
			     col.names = c("qseqid", "qlen", "sseqid", "slen", "length", "qstart", 
					   "qend", "sstart", "send", "bitscore", "evalue", "pident"))

# We want only the first alignment from a given query-subject combination
blastn_results <- blastn_results[!duplicated(blastn_results[, c("qseqid", "sseqid")]), ]

# Just a sanity check that all query sequences IDs have match within the VCF file
stopifnot(all(blastn_results$qseqid %in% rownames(vcf)))

# Executing the grep/awk script to generate metadata on the TE sequences
# DEPENDENCY : extract_te_metadata.sh
system("./extract_te_metadata.sh")

# Reading in the data that was generated
te_metadata <- read.table("te_metadata.txt", header = TRUE, stringsAsFactors = FALSE)

# --- This chunk filters the alignments found with blast by removing
# --- those for which either the query sequence length or subject
# --- sequence length is less than 80% of the total alignment length

# Let us add a column for the ratio of alignment length to query length
blastn_results$query_ratio <- blastn_results$length / blastn_results$qlen

# And another for the ratio of alignment length to subject length
blastn_results$subject_ratio <- blastn_results$length / blastn_results$slen

# We take the minimum of these two as the overlap column
blastn_results$overlap <- pmin(blastn_results$query_ratio, blastn_results$subject_ratio)

# Creating a vector of the names of query sequences for which this holds true
te_matches <- tapply(blastn_results$overlap, blastn_results$qseqid, function(x) any(x >= 0.8))
te_matches <- names(te_matches)[te_matches]

# Adding the data on SVTYPE and TE type to the blastn_results data.frame before filtering it
# Extracting the TE type from the sseqid string
blastn_results$te_type <- substr(blastn_results$sseqid, 6, 8)

# Classifying the queries as insertions or deletions. For this, we link to the VCF file
svtypes <- info(vcf)$SVTYPE
names(svtypes) <- rownames(vcf)
blastn_results$svtype <- svtypes[blastn_results$qseqid]

# Extracting a data.frame with the query names in te_matches and overlap >= 0.8
blastn_matches <- blastn_results[blastn_results$qseqid %in% te_matches & blastn_results$overlap >= 0.8, ]

# Filtering to keep only the best alignment for a given query sequence
blastn_matches <- blastn_matches[!duplicated(blastn_matches$qseqid), ]

# We can group them into meaningful groups for the comparisons that will be done downstream
te_lookup <- c("DHH" = "DNA", 
	       "DTA" = "DNA", 
	       "DTC" = "DNA", 
	       "DTH" = "DNA", 
	       "DTM" = "DNA", 
	       "DTT" = "DNA", 
	       "RIL" = "NON-LTR", 
	       "RIu" = "NON-LTR", 
	       "RLC" = "Copia", 
	       "RLG" = "Gypsy")

# Creating a new column with this classification
blastn_matches$te_group <- te_lookup[blastn_matches$te_type]

# --- This chunk adds some metadata to the blastn_matches data.frame
# --- It then also filters out non-polymorphic TE insertions and
# --- writes them to disk so they can be used for other analyses

# We can extract allele frequencies from the VCF
# Let us first extract only the variants in the vcf that correspond to known polymorphic TEs
vcf <- vcf[rownames(vcf) %in% te_matches, ]

# Let us do a sanity check
stopifnot(nrow(vcf) == nrow(blastn_matches))

# Let's extract the ALT allele frequency
alt_freq <- info(vcf)$AF
names(alt_freq) <- rownames(vcf)
# From it we can get the REF allele frequency
ref_freq <- 1 - alt_freq
# We can add those columns to the blastn_matches data.frame
blastn_matches$ref_freq <- as.numeric(ref_freq[blastn_matches$qseqid])
blastn_matches$alt_freq <- as.numeric(alt_freq[blastn_matches$qseqid])

# We can calculate the frequency of TE presence at the locus:
# If it is a deletion, then the TE frequency is the REF frequency
# If it is an insertion, then the TE frequency is the ALT frequency
blastn_matches$te_freq <- ifelse(blastn_matches$svtype == "DEL", blastn_matches$ref_freq, blastn_matches$alt_freq)

# We complete by adding the genomic positions of the polymorphic TEs to the blastn_matches data.frame
chrom <- as.character(seqnames(vcf))
names(chrom) <- rownames(vcf)

start <- as.numeric(start(vcf))
names(start) <- rownames(vcf)

end <- as.numeric(end(vcf))
names(end) <- rownames(vcf)

blastn_matches$chrom <- chrom[blastn_matches$qseqid]
blastn_matches$start <- start[blastn_matches$qseqid]
blastn_matches$end <- end[blastn_matches$qseqid]

# We extract the TEs that are polymorphic in the population (i.e. te_freq != 1 or 0)
polymorphic_tes <- blastn_matches[blastn_matches$te_freq != 1 & blastn_matches$te_freq != 0, ]

# Joining the metadata in the te_metadata data.frame with the one in polymorphic_tes
# First we need to generate the name column from the sseqid column
polymorphic_tes$name <- substring(polymorphic_tes$sseqid, 6)
stopifnot(all(polymorphic_tes$name %in% te_metadata$name))

# We reassign te_metadata to just the rows we need to cbind it to polymorphic_tes in the right order
stopifnot(!anyDuplicated(te_metadata$name))
rownames(te_metadata) <- te_metadata$name
te_metadata <- te_metadata[polymorphic_tes$name, ]
stopifnot(identical(polymorphic_tes$name, te_metadata$name))

# And now we cbind the columns together, except the name column from te_metadata because we already have it in polymorphic_tes
polymorphic_tes <- cbind(polymorphic_tes, te_metadata[, -1])

# We save those polymorphic TEs to file so we can use thse data for downstream analyses
# OUTPUT : te_analysis/polymorphic_tes.tsv
write.table(polymorphic_tes, file = "polymorphic_tes.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

