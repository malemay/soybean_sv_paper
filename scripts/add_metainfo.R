# File initially created on October 5, 2020

# Adding the header meta-information lines relative to the INFO tags added
# by the refine_breakpoints() function to a VCF.

# Will make both the VCF more robust and informative, and enable parsing
# with bioinformatic tools such as bcftools

# The function reads a vcf file, finds the last ##INFO line, and adds the meta-
# information after those lines

add_metainfo <- function(input_vcf, output_vcf) {
	# Reading the vcf lines
	vcf_lines <- scan(input_vcf, what = character(), sep = "\n", quiet = TRUE)

	# Finding the last ##INFO line
	last_info <- max(grep("^##INFO", vcf_lines))

	# Creating a character vector with the lines to add
	info_lines <- c('##INFO=<ID=MAX_SVLEN,Number=0,Type=Flag,Description="Breakpoints not refined because SVLEN exceeds the threshold set">',
			'##INFO=<ID=NO_ALIGNMENT,Number=0,Type=Flag,Description="Breakpoints not refined because no alignement was produced">',
			'##INFO=<ID=GATHER_DATA_FAILED,Number=0,Type=Flag,Description="Breakpoints not refined because the gathering alignment data step failed">',
			'##INFO=<ID=BELOW_THRESHOLDS,Number=0,Type=Flag,Description="Breakpoints not refined because quality thresholds were not met">',
			'##INFO=<ID=REALIGNED,Number=0,Type=Flag,Description="Breakpoints refined by the refine_breakpoints() R function">',
			'##INFO=<ID=ORIGINAL_START,Number=1,Type=Integer,Description="Start position of the variant before breakpoints were refined">',
			'##INFO=<ID=ORIGINAL_LEN,Number=1,Type=Integer,Description="Length of the variant before breakpoints were refined">',
			'##INFO=<ID=ORIGINAL_ALT,Number=1,Type=String,Description="ALT allele of the variant before breakpoints were refined (only applies for insertions)">',
			'##INFO=<ID=CONTIG_LEN,Number=1,Type=Integer,Description="Length of the contig used for breakpoint refinement">',
			'##INFO=<ID=IDENTITY,Number=1,Type=Integer,Description="Percent identity of the alignment used to refine the breakpoints">',
			'##INFO=<ID=GAPS,Number=1,Type=Integer,Description="Percentage of gaps in the alignment used to refine the breakpoints">',
			'##INFO=<ID=OVERLAP,Number=1,Type=Float,Description="Reciprocal overlap between the original and breakpoint-refined variant (only applies to deletions)">',
			'##INFO=<ID=DISTANCE_RATIO,Number=1,Type=Float,Description="Ratio of the Levenshtein distance between the original and breakpoint-refined allele to the length of the longest allele (only applies to insertions)">',
			'##INFO=<ID=OFFSET,Number=1,Type=Integer,Description="Absolute difference in the position of the insertion between the original allele and the breakpoint-refined allele, in nucleotides (only applies to insertions)">')

	# Updating the vcf_lines with those metainfo lines
	vcf_lines <- c(vcf_lines[1:last_info], info_lines, vcf_lines[(last_info + 1):length(vcf_lines)])

	# Writing the lines to the output file
	cat(vcf_lines, file = output_vcf, sep = "\n")

	# Returning NULL because the usable output is to file
	return(invisible(NULL))
}

