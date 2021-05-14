# File initially created on October 2, 2020

# This function takes the name of a vcf file as input and processes 
# the input line by line to write the output back to another vcf
# file.

# Many other functions work under the hood to allow this function
# to process the vcf lines in parallel and call the external commands
# that produce the alignment:

# - update_breakpoints()

refine_breakpoints <- function(input_vcf, output_vcf, ncores = 1,
			       reference_window, reads_window,
			       min_overlap, min_identity, max_gaps,
			       max_distance, max_offset, max_svlen,
			       age_script, samtools, minimap2,
			       age, wtdbg2, wtpoa_cns, refgenome,
			       nanopore_bam) {

	# Reading the input vcf_file using the scan function
	vcf_lines <- scan(input_vcf, what = character(), sep = "\n", quiet = TRUE)

	# Extracting the header lines and the body lines of the vcf
	header_lines <- grep("^#", vcf_lines, value = TRUE)
	body_lines <- grep("^#", vcf_lines, value = TRUE, invert = TRUE)

	# A sanity check to make sure that the splitting worked properly
	stopifnot(length(vcf_lines) == length(header_lines) + length(body_lines))

	# Here we process the header lines to add the tags that our function will add
	#
	# -
	# TODO
	#
	# -

	# We also write the header to the output file
	output_file <- file(output_vcf, open = "w")
	on.exit(close(output_file), add = TRUE)
	cat(header_lines, file = output_file, sep = "\n")

	# We coerce the lines of the vcf body to a list for processing with parallel::mclapply()
	body_lines <- as.list(body_lines)

	# Processing with mclapply()
	realigned_lines <- parallel::mclapply(body_lines, FUN = update_breakpoints, 
					      mc.cores = ncores,
					      reference_window = reference_window,
					      reads_window = reads_window,
					      min_overlap = min_overlap,
					      min_identity = min_identity,
					      max_gaps = max_gaps,
					      max_distance = max_distance,
					      max_offset = max_offset,
					      max_svlen = max_svlen,
					      age_script = age_script,
					      samtools = samtools,
					      minimap2 = minimap2,
					      age = age,
					      wtdbg2 = wtdbg2,
					      wtpoa_cns = wtpoa_cns,
					      refgenome = refgenome,
					      nanopore_bam = nanopore_bam)

	# Turning the list back into a character vector
	realigned_lines <- as.character(realigned_lines)

	# Writing the results to the output file
	cat(realigned_lines, file = output_file, sep = "\n")

	# No useful value is returned because the output is written to file
	return(invisible(NULL))
}

