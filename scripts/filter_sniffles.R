# File initially created on September 17, 2020

# This function was adapted from the function read_filter_vcf used as part of the sveval
# benchmarks. It uses sveval::readSVvcf's quick and convenient C++ parsing and uses this
# information to filter out variants. Filters applied to the GRanges object returned by
# readSVvcf are then used to selectively remove vcf lines parsed using R's scan function.

# WARNING: The argument min_sv_size is applied AFTER readSVvcf whereas the arugment
# min.sv.size is applied internally in the C++ code invoked by readSVvcf. In the context
# of this function, no filtering based on SV size is wanted in the C++ step because there
# has to be a direct correspondence between line numbers in the lines parsed by R and those
# parsed by readSVvcf.

filter_sniffles <- function(vcf.file, keep_hets = TRUE, max_N_del = NULL, keep_homref = FALSE,
			    chrom_pattern = NULL, remove_N_ins = FALSE, N_ins_range = NULL, refseq = NULL,
			    min_sv_size = 50, max_sv_size = 500000, output_file = NULL,
			    keep.ins.seq = FALSE, keep.ref.seq = FALSE, sample.name = "", 
			    qual.field = c("GQ", "QUAL"), other.field = NULL,
			    check.inv = FALSE, keep.ids = FALSE, nocalls = FALSE,
			    out.fmt = c("gr", "df", "vcf"), min.sv.size = 10) {

	# Checking that an output file has been provided and opening it for writing
	if(is.null(output_file)) stop("output_file must be provided")
	output <- file(output_file, open = "w")
	on.exit(close(output))

	# Reading the lines from the sniffles file with R
        message("Reading ", vcf.file, " using R's scan() function")
	r_lines <- scan(vcf.file, what = character(), sep = "\n")

	# Outputting the header lines to file and keeping only the data lines in memory
	# Getting the indices of header lines
	header_lines <- grep("^#", r_lines)
	# Checking that header lines are located as expected
	if(! 1 %in% header_lines || ! all(diff(header_lines) == 1)) stop("Unexpected header format")
	# Outputting the header lines to file
	cat(r_lines[header_lines], sep = "\n", file = output)
	# Removing the header_lines from r_lines
	r_lines <- r_lines[-header_lines]

	# Breakends are removed from the R lines
	breakends <- grepl("SVTYPE=BND", r_lines)
	message("Removing ", sum(breakends), " SVs with SVTYPE=BND")
	r_lines <- r_lines[!breakends]
	
	# Reading the vcf file with sveval:::read_vcf_cpp with the specified parameters
	message("Reading ", vcf.file, " using sveval::read_vcf_cpp")
	svs <- sveval:::read_vcf_cpp(vcf.file, use_gz = FALSE, sample_name = sample.name,
				     shorten_ref = !keep.ref.seq, shorten_alt = !keep.ins.seq,
				     check_inv = check.inv, gq_field = qual.field[1],
				     keep_nocalls = nocalls, other_field = other.field,
				     min_sv_size = min.sv.size)
	
	# Checking that the number of variants is the same in both files
	if(length(r_lines) != nrow(svs)) stop("r_lines and svs do not have the same length")

	# Also checking that the chromomsome and starting positions are the same
	r_chr <- sapply(strsplit(r_lines, "\t"), function(x) x[1])
	r_start <- sapply(strsplit(r_lines, "\t"), function(x) x[2])

	if(! all(r_chr == svs$seqnames) || ! all(r_start == svs$start)) {
		stop("Data differ between scan()ned data and data read with read_vcf_cpp")
	}

	# Log info
	message("Initial SV set contains ", nrow(svs), " SVs")

	# Creating a column with line numbers so the R file can be output at the end
	svs$line_number <- 1:nrow(svs)

	# Removing variants that have a negative width ; otherwise no conversion to GRanges possible
	if(any(neg_width <- svs$start > svs$end)) {
		message("From readSVvcf: ", sum(neg_width), " variants removed due to start position > end position")
		svs <- svs[!neg_width, ]
	}

	# Converting do GRanges
	svs <- GenomicRanges::makeGRangesFromDataFrame(svs, keep.extra.columns=TRUE)

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

	# Filtering out insertions that have N within N_ins_range of their location
	if(!is.null(N_ins_range)) {
		# Checking that the reference is provided and exists
		stopifnot(!is.null(refseq) && file.exists(refseq))

		# Extracting the maximum position to extract for each SV depending on its chromosome
		chr_ranges <- Rsamtools::scanFaIndex(refseq)
		max_coord <- end(chr_ranges)[as.numeric(match(seqnames(svs), seqnames(chr_ranges)))]
		max_coord <- as.numeric(max_coord)

		# Extracting the flanking nucleotides of all variants
		flank_seq <- Rsamtools::scanFa(refseq, GRanges(seqnames = seqnames(svs),
							       IRanges(start = pmax(start(svs) - N_ins_range, 1),
								       end = pmin(start(svs) + N_ins_range, max_coord))))

		# Removing insertions containing Ns in flanking seq
		flanking_N <- svs$type == "INS" & grepl("N", as.character(flank_seq), fixed = TRUE)
		message(sum(flanking_N), " insertions removed due to N in flanking ", N_ins_range, " nucleotides")
		svs <- svs[!flanking_N]
	}

	# Remove SVs smaller than min_sv_size
	if(!is.null(min_sv_size)) {
		min_sv <- svs$size < min_sv_size
		message(sum(min_sv), " variants smaller than ", min_sv_size, " nucleotides removed.")
		svs <- svs[!min_sv]
	}

	# Remove SVs larger than max_sv_size
	if(!is.null(max_sv_size)) {
		max_sv <- svs$size > max_sv_size
		message(sum(max_sv), " variants larger than ", max_sv_size, " nucleotides removed.")
		svs <- svs[!max_sv]
	}

	# Remove insertions larger than 46340 to avoid integer overflow during insertion comparison
	if(any(long_ins <- svs$type == "INS" & nchar(svs$alt) > 46340)) {
		message(sum(long_ins), " insertion alleles larger than 46340 nucleotides removed")
		svs <- svs[!long_ins]
	}

	# Filterting deletions/duplications/inversions whose sequence contains a proportion of > N_max N nucleotides
	if(!is.null(max_N_del)) {

		# Changing the <DEL> reference sequences to the real sequence
		if(any(svs$alt %in% c("<DEL>", "<INV>", "<DUP>"))) {
			# Checking that refseq is provide and exists
			stopifnot(!is.null(refseq) && file.exists(refseq))

			# Removing the deletions/inversions/duplications that are larger than 500,000 nucleotides (if any not removed yet)
			large_del <- svs$size > 500000 & svs$type %in% c("DEL", "INV", "DUP")

			if(any(large_del)) {
				message(sum(large_del), " variants larger than 500,000 nucleotides removed")
				svs <- svs[!large_del]
			}

			# Getting a vector of indices of the <DEL> alleles
			del_alleles <- which(svs$alt %in% c("<DEL>", "<INV>", "<DUP>")) 
			message("Found ", length(del_alleles), " variants to query in ", refseq)
			svs$ref[del_alleles] <- Rsamtools::scanFa(refseq, svs[del_alleles])
		}

		N_deletions <- svs$type %in% c("DEL", "INV", "DUP") & stringr::str_count(svs$ref, "N")/nchar(svs$ref) > max_N_del
		message(sum(N_deletions & svs$type == "DEL"), " deletions removed due to a proportion of N > ", max_N_del)
		message(sum(N_deletions & svs$type == "INV"), " inversions removed due to a proportion of N > ", max_N_del)
		message(sum(N_deletions & svs$type == "DUP"), " duplications removed due to a proportion of N > ", max_N_del)
		svs <- svs[!N_deletions]

		# Replacing the <DEL> deletions back to N in order to spare memory/disk space
		if(any(svs$alt %in% c("<DEL>", "<INV>", "<DUP>"))) svs$ref[svs$alt %in% c("<DEL>", "<INV>", "<DUP>")] <- "N"	
	}

	# Selecting the lines to write
	r_lines <- r_lines[svs$line_number]

	# Last filtering step be removing "UNRESOLVED" variants
	unresolved <- grepl("UNRESOLVED", r_lines)
	message(sum(unresolved), " UNRESOLVED variants removed")
	r_lines <- r_lines[!unresolved]

	# Log info
	message("Final SV set contains ", length(r_lines), " SVs")

	# Writing final dataset to file
	message("Writing filtered vcf file to ", output_file)
	cat(r_lines, sep = "\n", file = output)

	return(invisible(svs))
}

