# File initially created on September 10, 2020

# This function reads a vcf file using sveval::readSVvcf (arguments can be passed to it)
# and then optionally applies some filters to customize the structural variants that
# will be compared by svevalOl

read_filter_vcf  <-  function(vcf.file, keep_hets = TRUE, max_N_del = NULL, keep_homref = FALSE,
			      chrom_pattern = NULL, remove_N_ins = FALSE, refseq = NULL,
			      keep.ins.seq = FALSE, keep.ref.seq = FALSE, sample.name = "", 
			      qual.field = c("GQ", "QUAL"), other.field = NULL,
			      check.inv = FALSE, keep.ids = FALSE, nocalls = FALSE,
			      out.fmt = c("gr", "df", "vcf"), min.sv.size = 10) {

	# Log info
	message("Reading sample ", sample.name, " from vcf file ", vcf.file)

	# Reading the vcf file with the specified parameters
	svs  <- readSVvcf(vcf.file, keep.ins.seq = keep.ins.seq,
			  keep.ref.seq = keep.ref.seq, sample.name = sample.name,
			  qual.field = qual.field, other.field = other.field,
			  check.inv = check.inv, keep.ids = keep.ids,
			  nocalls = nocalls, out.fmt = out.fmt,
			  min.sv.size = min.sv.size)
	
	# Log info
	message("Initial SV set contains ", length(svs), " SVs")

	# Removing SVs that are not insertions or deletions
	other_types <- setdiff(unique(svs$type), c("INS", "DEL", "INV", "DUP"))

	# Looping to remove these types one by one
	if(length(other_types)) {
		for(i in other_types) {
			to_remove <- svs$type == i
			message(sum(to_remove), " variants of type ", i, " removed")
			svs <- svs[!to_remove]
		}
	}

	# Removing reference allele calls
	if(!keep_homref) {
		# Getting a vector of homozygous reference calls
		homref <- svs$ac == 0
		# Printing the number of matches removed
		message(sum(homref), " homozygous reference calls removed")
		# Filtering out
		svs <- svs[!homref]
	}

	# Removing heterozygous calls
	if(!keep_hets) {
		# Getting a vector of heterozygous calls
		het <- svs$ac == 1
		# Printing the number of matches removed
		message(sum(het), " heterozygous calls removed")
		# Filtering out
		svs <- svs[!het]
	}

	# Removing the "*" alternative alleles because they do not represent true deletions
	star_alleles <- svs$alt == "*"
	message(sum(star_alleles), ' "*" alleles removed')
	svs <- svs[!star_alleles]

	# Filtering out the seqnames not matching chrom_pattern if not NULL
	if(!is.null(chrom_pattern)) {
		# Getting a logical vector of the positions matching pattern
		chrom_match <- grepl(chrom_pattern, GenomicRanges::seqnames(svs))
		# Printing the number of matches removed
		message(sum(!chrom_match), " alleles not matching ", chrom_pattern, " removed")
		# Filtering them out
		svs <- svs[chrom_match]
	}

	# Filtering out the insertions that have "N" in their sequence
	if(remove_N_ins) {
		# Removing insertions with an undefined sequence
		imprecise_ins <- svs$alt == "<INS>"
		message(sum(imprecise_ins), " insertions removed due to unknown sequence (<INS>)")
		svs <- svs[!imprecise_ins]

		# Removing insertions with any N in the ALT sequence
		N_insertions <- svs$type == "INS" & grepl("N", svs$alt)
		message(sum(N_insertions), " insertions removed due to ALT containing at least one N")
		svs <- svs[!N_insertions]

		# Removing insertions with > 1 N in the reference sequence
		N_insertions <- svs$type == "INS" & stringr::str_count(svs$ref, "N") > 1
		message(sum(N_insertions), " removed due to REF containing > 1 N")
	}

	# Remove insertions larger than 46340 to avoid integer overflow during insertion comparison
	if(any(long_ins <- svs$type == "INS" & nchar(svs$alt) > 46340)) {
		message(sum(long_ins), " insertion alleles larger than 46340 nucleotides removed")
		svs <- svs[!long_ins]
	}

	# Filterting deletions whose sequence contains a proportion of > N_max N nucleotides
	if(!is.null(max_N_del)) {

		# Changing the <DEL> reference sequences to the real sequence
		if(any(svs$alt == "<DEL>")) {
			# Checking that refseq is provided and exists
			stopifnot(!is.null(refseq) && file.exists(refseq))

			# Removing the deletions that are larger than 500,000 nucleotides
			large_del <- svs$size > 500000 & svs$type == "DEL"
			message(sum(large_del), " deletions larger than 500,000 nucleotides removed")
			svs <- svs[!large_del]

			# Getting a vector of indices of the <DEL> alleles
			del_alleles <- which(svs$alt == "<DEL>")
			message("Found ", length(del_alleles), " deletions to query in ", refseq)
			svs$ref[del_alleles] <- Rsamtools::scanFa(refseq, svs[del_alleles])
		}

		N_deletions <- svs$type == "DEL" & stringr::str_count(svs$ref, "N")/nchar(svs$ref) > max_N_del
		message(sum(N_deletions), " deletions removed due to a proportion of N > ", max_N_del)
		svs <- svs[!N_deletions]

		# Replacing the <DEL> deletions back to N in order to spare memory/disk space
		if(any(svs$alt == "<DEL>")) svs$ref[svs$alt == "<DEL>"] <- "N"	
	}

	# Log info
	message("Final SV set contains ", length(svs), " SVs")

	return(svs)
}

