#!/prg/R/4.0/bin/Rscript

# File initially created on April 21, 2021

# This code performs the Gene Ontology analysis to find terms
# that are over- or underrepresented among genes that are overlapped
# by high-frequency (> 0.75 SV frequency) SVs.

# The analysis will be carried out using the GOstats package
# using the conditional test. A Bonferroni correction will be
# applied based on the number of terms tested to test for
# significance

# Loading the required packages
library(GenomicRanges)
library(GOstats)
library(GSEABase)
library(PFAM.db)

# Loading the GRanges object of genes
# DEPENDENCY : genes.RData
load("genes.RData")

#----- Preparing the list of genes that overlap high-frequency SVs
# DEPENDENCY : overlap_data.txt
overlap_data <- read.table("overlap_data.txt", header = TRUE, stringsAsFactors = FALSE)

# Extracting the SVs that overlap genes (either coding or non-coding sequence)
frequent_gene_svs <- overlap_data[overlap_data$overlap == "cds", ]

# Then only extracting the ones with an allele frequency >= 0.5
frequent_gene_svs <- frequent_gene_svs[frequent_gene_svs$af >= 0.5, ]

# And then removing the ones that are fixed in the population because they might represent errors in the reference
frequent_gene_svs <- frequent_gene_svs[frequent_gene_svs$af != 1, ]

# Coercing to a GRanges so we can check the overlap with genes
frequent_gene_svs <- makeGRangesFromDataFrame(frequent_gene_svs)

# Let us extract the list of genes that overlap these
sv_genes <- subsetByOverlaps(genes, frequent_gene_svs)
sv_genes <- sv_genes$gene_id
sv_genes <- toupper(sv_genes)

# Reading the annotation in the file retrieved from Soybease on April 20, 2021
# First using a system call to simplify the file to only what we need
# DEPENDENCY : soybase_genome_annotation_v4.0_04-20-2021.txt
system("tail -n +14 soybase_genome_annotation_v4.0_04-20-2021.txt | cut -f 1,12-18 > gmax_go_annotations.txt")
# annotations <- read.table("gmax_go_annotations.txt", header = TRUE, sep = "\t")

# read.table does not seem able to read the data properly, so I will have to set the data.frame manually
fields <- strsplit(scan("gmax_go_annotations.txt", what = character(), sep = "\n", quiet = TRUE), "\t")

# Checking that the number of fields is the same for all lines
table(lengths(fields))
# 
#     8 
# 52873 
# all good

# Now coercing this data into a data.frame and doing some formatting
gmax_annotations <- as.data.frame(do.call("rbind", fields), stringsAsFactors = FALSE)
colnames(gmax_annotations) <- gmax_annotations[1, ]
gmax_annotations <- gmax_annotations[-1, ]
names(gmax_annotations) <- c("gene", "go_bp_id", "go_bp_desc", "go_mf_id", "go_mf_desc", "go_cc_id", "go_cc_desc", "pfam_id")

# Transforming all the "NA" strings to real NA values
for(i in 1:ncol(gmax_annotations)) gmax_annotations[[i]][gmax_annotations[[i]] == "NA"] <- NA

# Formatting the gene_names in all caps
gmax_annotations$gene <- toupper(gmax_annotations$gene)

# We only keep genes that are located on chromosomes (i.e. not on unanchored scaffolds)
# We can select them grom the ID
gmax_annotations <- gmax_annotations[grepl("GLYMA\\.[0-9]{2}G[0-9]{6}", gmax_annotations$gene), ]

# This should correspond to the number of genes in the GRanges object loaded above
stopifnot(nrow(gmax_annotations) == length(genes))

# Checking the overlap between our sv_genes vector and the annotations
stopifnot(all(sv_genes %in% gmax_annotations$gene))

# Creating a function that takes a vector of gene IDs and a vector of GO IDs and generates
# a GOAllFrame object ready for input to GeneSetCollection
annotation_to_goframe <- function(gene_ids, go_ids) {

	# Removing the entries for which there is no GO annotation	
	gene_ids <- gene_ids[!is.na(go_ids)]
	go_ids   <- go_ids[!is.na(go_ids)]

	# Splitting the go_ids using the space between them
	go_ids <- strsplit(go_ids, " ")

	# Using the number of GO terms per gene to generate a proper vector for genes
	gene_ids <- rep(gene_ids, times = lengths(go_ids))

	# Now we can generate the final data.frame
	output <- data.frame(go_id = unlist(go_ids),
			     evidence = rep("IEA", length(gene_ids)),
			     gene_id = gene_ids,
			     stringsAsFactors = FALSE)

	GOAllFrame(GOFrame(output, organism = "Glycine max"))
}

# I am more interested in the Biological Process ontology so I will only test this one
bp_goframe <- annotation_to_goframe(gmax_annotations$gene, gmax_annotations$go_bp_id)

# We need to generate a geneSetCollection object
bp_gsc <- GeneSetCollection(bp_goframe, setType = GOCollection())

# I keep only sv_genes for which there is a BP annotation
bp_genes <- sv_genes[sv_genes %in% gmax_annotations[!is.na(gmax_annotations$go_bp_id), "gene"]]

bp_over_params <- GSEAGOHyperGParams(name="Biological Process overexpression parameters",
				     geneSetCollection=bp_gsc,
				     geneIds = bp_genes,
				     universeGeneIds = gmax_annotations[!is.na(gmax_annotations$go_bp_id), "gene"],
				     ontology = "BP",
				     pvalueCutoff = 0.05,
				     conditional = TRUE,
				     testDirection = "over")

# Computing the test using those parameters
bp_over_test <- hyperGTest(bp_over_params)
bp_over_summary <- summary(bp_over_test, categorySize = 20, pvalue = 1)
bp_over_summary$Pvalue <- bp_over_summary$Pvalue * nrow(bp_over_summary)

# I also want to test for GO term depletion
bp_under_params <- GSEAGOHyperGParams(name="Biological Process overexpression parameters",
				      geneSetCollection=bp_gsc,
				      geneIds = bp_genes,
				      universeGeneIds = gmax_annotations[!is.na(gmax_annotations$go_bp_id), "gene"],
				      ontology = "BP",
				      pvalueCutoff = 0.05,
				      conditional = TRUE,
				      testDirection = "under")

# Computing the test using those parameters
bp_under_test <- hyperGTest(bp_under_params)
bp_under_summary <- summary(bp_under_test, categorySize = 20, pvalue = 1)
bp_under_summary$Pvalue <- bp_under_summary$Pvalue * nrow(bp_under_summary)

# The PFAM terms analysis is not supported for Glycine max because I do not have
# an annotation package and there is no function from GSEABase that allows me to
# generate suitable input. Therefore, I will implement a homemade hypergeometric 
# test to test for PFAM terms enrichment. 

# There appears to be an internal function in the Category package (Category:::.doHyperGInternal)
# that I could take advantage of to create my own testing function

# Getting a PFAM data.frame using a refactored function similar to the one used for creating a GOFrame object
annotation_to_pfamframe <- function(gene_ids, pfam_ids) {

	# Removing the entries for which there is no PFAM annotation	
	gene_ids <- gene_ids[!is.na(pfam_ids)]
	pfam_ids   <- pfam_ids[!is.na(pfam_ids)]

	# Splitting the pfam_ids using the space between them
	pfam_ids <- strsplit(pfam_ids, " ")

	# Using the number of GO terms per gene to generate a proper vector for genes
	gene_ids <- rep(gene_ids, times = lengths(pfam_ids))

	# Now we can generate the final data.frame
	output <- data.frame(pfam_id = unlist(pfam_ids),
			     gene_id = gene_ids,
			     stringsAsFactors = FALSE)

	output
}

# Applying the function to generate the input data for the test
pfam_frame <- annotation_to_pfamframe(gmax_annotations$gene, gmax_annotations$pfam_id)

# Some domains are present more than once for a given gene, so we remove duplicates
pfam_frame <- pfam_frame[!duplicated(pfam_frame), ]

# Testing if we can retrieve all the PFAM IDs from the PFAM.db annotation package
pfam_ids <- unique(pfam_frame$pfam_id)
pfam_domains <- unlist(as.list(PFAM.db::PFAMDE)[pfam_ids])
stopifnot(all(!is.na(pfam_domains)))

# Creating a function that tests enrichment similarly to what hyperGTest does
# gene_ids : the genes that were found by the analysis
# pfam_data: A data.frame with a pfam_id column and a gene_id column
#            that describes the PFAM annotation for the gene universe considered
# categorySize: the minimum number of genes in the universe corresponding to this
#               PFAM domain for it to be considered
# over: whether we are to test for overrepresentation (TRUE) or underrepresentation (FALSE)
pfam_hyperGTest <- function(gene_ids, pfam_data, categorySize = 10, over = TRUE) {

	# Extracting the part of the data.frame corresponding to the gene that were selected
	gene_data <- pfam_data[pfam_data$gene_id %in% gene_ids, ]
	# Counting the number of instances of each PFAM ID within this set of genes (observed counts)
	pfam_table <- table(gene_data$pfam_id)
	# Creating a data.frame out of these data and adding a column with the count of each PFAM ID in the gene universe
	pfam_to_test <- data.frame(id = names(pfam_table),
				   observed = as.numeric(pfam_table),
				   in_universe = as.numeric(table(pfam_data$pfam_id)[names(pfam_table)]),
				   stringsAsFactors = FALSE)

	# If testing for under-representation, we also need to add PFAM IDs for which the count in gene_ids is 0
	if(!over) {
		universe_table <- table(pfam_data$pfam_id)
		universe_table <- universe_table[!names(universe_table) %in% pfam_to_test$id]
		universe_df    <- data.frame(id = names(universe_table),
					     observed = 0,
					     in_universe = as.numeric(universe_table),
					     stringsAsFactors = FALSE)
		pfam_to_test <- rbind(pfam_to_test, universe_df)
	}

	# Removing the PFAM entries that have less than categorySize matches in the reference
	pfam_to_test <- pfam_to_test[pfam_to_test$in_universe >= categorySize, ]

	# Extracting the variables that will be used for all hypergeometric tests
	n_drawn    <- length(gene_ids)
	n_universe <- length(unique(pfam_data$gene_id))

	# Initializing an output data.frame
	output <- data.frame(PFAMID = character(),
			     Pvalue = numeric(),
			     OddsRatio = numeric(),
			     ExpCount = numeric(),
			     Count = numeric(),
			     Size = numeric(),
			     stringsAsFactors = FALSE)

	# Looping over all the PFAM IDs that need to be tested
	for(i in 1:nrow(pfam_to_test)) {
		# Extracting some values and computing the test
		n_observed <- pfam_to_test[i, "observed"]
		i_universe <- pfam_to_test[i, "in_universe"]
		i_results <- Category:::.doHyperGInternal(numW = i_universe, numB = n_universe - i_universe,
							  numDrawn = n_drawn, numWdrawn = n_observed,
							  over = over)

		# Adding the data for this PFAM entry to the output data.frame
		i_results <- data.frame(PFAMID = pfam_to_test[i, "id"],
					Pvalue = i_results$p,
					OddsRatio = i_results$odds,
					ExpCount = i_results$expected,
					Count = n_observed,
					Size = i_universe,
					stringsAsFactors = FALSE)

		output <- rbind(output, i_results)
	}

	# Adding the PFAM terms and reordering the output data.frame
	output$Term <- NA
	PFAM_terms <- as.list(PFAM.db::PFAMDE)
	
	for(i in 1:nrow(output)) {
		if(output[i, "PFAMID"] %in% names(PFAM_terms)) output[i, "Term"] <- PFAM_terms[[output[i, "PFAMID"]]]
	}

	output[order(output$Pvalue), ]
}

# Now we can test for over- and underrepresentation of PFAM terms using this function
# of PFAM terms within our dataset of coding sequences impacted by SVs
pfam_genes <- sv_genes[sv_genes %in% gmax_annotations[!is.na(gmax_annotations$pfam_id), "gene"]]

# Computing the test for over-representation
pfam_over_test <- pfam_hyperGTest(pfam_genes, pfam_frame, categorySize = 10, over = TRUE)
pfam_over_test$Pvalue <- pfam_over_test$Pvalue * nrow(pfam_over_test)
# Renaming the result for consistency with the naming scheme for the GO tests
pfam_over_summary <- pfam_over_test

# Computing the test for under-representation
pfam_under_test <- pfam_hyperGTest(pfam_genes, pfam_frame, categorySize = 10, over = FALSE)
pfam_under_test$Pvalue <- pfam_under_test$Pvalue * nrow(pfam_under_test)
# We should not get any significant hits for this test but we'll test jsut in case in changes over time
stopifnot(nrow(pfam_under_test) != 0)

# Saving the objects from the GO analysis to file
save(bp_over_summary, file = "bp_over_summary.RData")
save(bp_under_summary, file = "bp_under_summary.RData")
save(pfam_over_summary, file = "pfam_over_summary.RData")

