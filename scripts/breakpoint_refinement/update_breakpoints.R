# File initially created on Friday, October 2, 2020

# This function accepts a vcf line of a structural variant as input
# and returns a line with the vcf coordinates optionally updated following
# realignment with AGE.

# The function is meant to be called via the refine_breakpoints() function
# which will apply it in parallel to every line of a vcf file using mclapply()

# This function requires the following functions to be available in the globalenv()
# - parse_svinfo() to parse the informations from the vcf line
# - call_age() which is a wrapper around the bash script age_realign.sh
# - parse_age() which gets in the information from the age file into R
# - gather_align_data() which makes a data.frame from the age and svinfo data

update_breakpoints <- function(vcf_line, reference_window, reads_window,
			       min_overlap, min_identity, max_gaps,
			       max_distance, max_offset, max_svlen,
			       age_script, samtools, minimap2, age,
			       wtdbg2, wtpoa_cns, refgenome, nanopore_bam) {

	# Parsing the information from the vcf line
	svinfo <- parse_svinfo(vcf_line)

	# We do not realign duplications and inversions; these are returned unmodified
	if(svinfo$svtype %in%  c("DUP", "INV")) {
		return(vcf_line)
	}

	# At this point there should only be insertions and deletions left
	stopifnot(svinfo$svtype %in% c("DEL", "INS"))

	# We do not process this SV if if is larger than max_svlen
	if(svinfo$svlen > max_svlen) {
		warning("Following line not processed, SV length > ", max_svlen, "\n", vcf_line)
		sline <- strsplit(vcf_line, "\\s+")[[1]]
		sline[8] <- paste0(sline[8], ";MAX_SVLEN")
		return(paste0(sline, collapse = "\t"))
	}

	# Creating two boolean values to make the code clearer
	del <- svinfo$svtype == "DEL"
	ins <- svinfo$svtype == "INS"

	# Passing the information to the call_age() function
	# This function returns a character vector of length 2 with
	# the information needed for the rest of the function :
	# - [1] the name of the .age.txt file created
	# - [2] the name of the .ontpolish.fa file containing the assembled and polished contig
	age_files <- call_age(chr = svinfo$chr, svtype = svinfo$svtype, start = svinfo$start,
			      svlen = svinfo$svlen, reference_window = reference_window, 
			      reads_window = reads_window, age_script = age_script,
			      samtools = samtools, minimap2 = minimap2, age = age,
			      wtdbg2 = wtdbg2, wtpoa_cns = wtpoa_cns, refgenome = refgenome,
			      nanopore_bam = nanopore_bam)

	# Parsing the age file
	age_data <- parse_age(age_files[1])

	# We check if age returned anything (will be NULL in that case)
	if(is.null(age_data)) {
		# In that case we return the vcf line (almost) unmodified and remove files
		file.remove(age_files[1], age_files[2])
		sline <- strsplit(vcf_line, "\\s+")[[1]]
		sline[8] <- paste0(sline[8], ";NO_ALIGNMENT")
		return(paste0(sline, collapse = "\t"))
	}

	# Otherwise we generate a data.frame with the data pertaining to that line
	align_data <- gather_align_data(svinfo, age_data, reference_window, age_files[2])

	# We no longer need the age and contig files
	file.remove(age_files[1], age_files[2])

	# We also remove the .fai index if it exists
	if(file.exists(fai <- paste0(age_files[2], ".fai"))) file.remove(fai)

	# We check if the call to gather_align_data returned anything useful
	if(is.null(align_data)) {
		sline <- strsplit(vcf_line, "\\s+")[[1]]
		sline[8] <- paste0(sline[8], ";GATHER_DATA_FAILED")
		return(paste0(sline, collapse = "\t"))
	}

	# Based on these data we decide how to realign the variant
	# There may be more than one row in the data.frame so we have to loop over
	# all of them ; we stop looping once we found a suitable one

	# This variable determines whether there has been realignment
	realigned <- FALSE

	for(i in 1:nrow(align_data)) {

		# Extract some data which is the same no matter if INS or DEL
		identity <- align_data[i, "identity"]
		gaps <- align_data[i, "gaps"]

		# We apply different filters if they are deletions or insertions
		if(del) {
			# We extract the reciprocal overlap
			rol <- align_data[i, "r_overlap"]

			if(rol >= min_overlap && identity >= min_identity && gaps <= max_gaps) {
				realigned <- TRUE
				irow <- i
				break
			}

		} else if (ins) {
			# Extracting the Levenshtein distance ratio and absolute offset
			dist_ratio <- align_data[i, "dist_ratio"]
			offset <- abs(align_data[i, "offset"])

			if(dist_ratio <= max_distance && offset <= max_offset && identity >= min_identity && gaps <= max_gaps) {
				realigned <- TRUE
				irow <- i
				break
			}
		}
	}

	# We return the original line if it is not to be updated
	if(!realigned) {
		sline <- strsplit(vcf_line, "\\s+")[[1]]
		sline[8] <- paste0(sline[8], ";BELOW_THRESHOLDS")
		return(paste0(sline, collapse = "\t"))
	}

	# Otherwise we need to update the vcf line
	# We first split the vcf line into its fields
	sline <- strsplit(vcf_line, "\\s+")[[1]]

	# We first update the start position, SVLEN, and ALT allele (if an insertion)
	sline[2] <- align_data[irow, "age_start"]
	if(ins) sline[5] <- align_data[irow, "age_alt"]

	# The SVLEN modification will depend on whether it is a deletion or an insertion
	if(del) {
		sline[8] <- sub("SVLEN=-\\d+",
				paste0("SVLEN=", as.character(-align_data[irow, "excise_ref_len"])),
				sline[8])
	} else if(ins) {
		sline[8] <- sub("SVLEN=\\d+",
				paste0("SVLEN=", as.character(align_data[irow, "excise_contig_len"])),
				sline[8])
	}

	# We also need to update INFO/END, which depends on whether it is INS or DEL
	# The way END is computed for deletions is not how I think about it
	# (I would have done start + length - 2)
	# But it is how Sniffles does, so I do the same out of consistency
	if(del) {
		sline[8] <- sub("([^[:alpha:]])END=\\d+",
				paste0("\\1END=", as.character(align_data[irow, "age_start"] + 
							    align_data[irow, "excise_ref_len"])),
				sline[8])
	} else if(ins) {
		sline[8] <- sub("([^[:alpha:]])END=\\d+",
				paste0("\\1END=", as.character(align_data[irow, "age_start"])),
				sline[8])
	}

	# Now we can add tags to the INFO field
	# The INFO field is the 8th one so we must update index 8
	sline[8] <- paste0(sline[8], ";REALIGNED")
	sline[8] <- paste0(sline[8], ";ORIGINAL_START=", svinfo$start)
	sline[8] <- paste0(sline[8], ";ORIGINAL_LEN=", ifelse(del, -svinfo$svlen, svinfo$svlen))
	if(ins) sline[8] <- paste0(sline[8], ";ORIGINAL_ALT=", svinfo$alt)
	sline[8] <- paste0(sline[8], ";CONTIG_LEN=", align_data[irow, "age_contiglen"])
	sline[8] <- paste0(sline[8], ";IDENTITY=", as.character(identity))
	sline[8] <- paste0(sline[8], ";GAPS=", as.character(gaps))
	if(del) sline[8] <- paste0(sline[8], ";OVERLAP=", as.character(rol))
	if(ins) sline[8] <- paste0(sline[8], ";DISTANCE_RATIO=", as.character(dist_ratio))
	if(ins) sline[8] <- paste0(sline[8], ";OFFSET=", as.character(offset))

	# By now we are able to return the re-assembled line
	return(paste0(sline, collapse = "\t"))
}

