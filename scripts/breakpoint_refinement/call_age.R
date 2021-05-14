# File initially created on Friday, October 2, 2020

# This function accepts arguments from update_breakpoints() in order
# to launch the SV realignment pipeline which involves the following steps:

# - Assembly of the Nanopore reads from the region of interest with wtdbg2
# - Polishing of the resulting assembly using the Nanopore reads themselvs
#   aligned to the assembled contig with minimap2
# - Alignement of the resulting contig to the reference genome using AGE

# This function is merely a wrapper around the bash script age_realign.sh

# The function creates some temporary files that are used by the bash script.
# It also manages those files by removing those that are no longer used at the
# end and returning the names of the alignment file (.age.txt) and contig fasta
# file (.ontpolish.fa) to the function that calls it so it knows where to find them

# This function also creates the command line that will be parsed using getopts
# by the bash script.

call_age <- function(chr, svtype, start, svlen, reference_window, reads_window,
		     age_script, samtools, minimap2, age, wtdbg2, wtpoa_cns,
		     refgenome, nanopore_bam) {

	# Creating temporary files that will be used by the bash script

	# The file to which the Nanopore reads will be written
	reads_fasta <- tempfile("reads_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(reads_fasta)
	on.exit(file.remove(reads_fasta), add = TRUE)

	# The file to which the assembled contig will be written
	contig_fasta <- tempfile("contig_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(contig_fasta)
	on.exit(file.remove(contig_fasta), add = TRUE)

	# The file to which the polished contig will be written
	polished_fasta <- tempfile("polished_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(polished_fasta)

	# The file to which the reference sequence in the region will be written
	reference_fasta <- tempfile("reference_fasta", tmpdir = getwd(), fileext = ".fa")
	file.create(reference_fasta)
	on.exit(file.remove(reference_fasta), add = TRUE)

	# The file to which the age alignments will be written
	age_file <- tempfile("age", tmpdir = getwd(), fileext = ".txt")
	file.create(age_file)

	# The prefix that fill be used
	prefix <- tempfile("prefix", tmpdir = "")
	prefix <- sub("/", "", prefix)

	command <- paste0(age_script,
			  " -s ", samtools,
			  " -m ", minimap2,
			  " -a ", age,
			  " -w ", wtdbg2,
			  " -c ", wtpoa_cns,
			  " -t ", svtype,
			  " -h ", chr,
			  " -e ", start,
			  " -l ", svlen,
			  " -r ", refgenome,
			  " -b ", nanopore_bam,
			  " -n ", reads_fasta,
			  " -f ", contig_fasta,
			  " -p ", polished_fasta,
			  " -z ", reference_fasta,
			  " -x ", age_file,
			  " -y ", prefix,
			  " -u ", reads_window,
			  " -d ", reference_window)

	# Launching the command itself
	system(command)

	# Returning the names of the files to be used by update_breakpoints()
	return(c(age_file, polished_fasta))
}

