# File initially created on Tuesday, September 29, 2020

# A function used to parse the output of the AGE aligner

# It takes the name of the file to be parsed (a string) as input and returns
# a list with the following elements:

# -
# -
# -
# -
# -

parse_age <- function(filename) {
  # Reading the whole file in memory using scan (file is overall small)
  age <- scan(filename, what = character(), sep = "\n")

  # Counting the number of alignments in the file
  # This is done by counting the number of matches to the regexp "^MATCH"
  block_start <- grep("^MATCH", age)
  n_align <- length(block_start)

  # Checking that any line matched
  if(n_align == 0) {
	  warning("No alignment returned for ", filename, ", returning NULL")
	  return(NULL)
  }


  # Otherwise checking the number of alignments
  # And selecting the first one if more than one is found
  if(n_align > 1) {

	  warning(n_align, " contigs aligned for ", filename, ", using only first one")

	  # Extracting the first part of the alignment (region between the two first ^MATCH)
	  age <- age[block_start[1]:(block_start[2] - 1)]
  }

  # Checking that there has been an excised fragment
  if(!any(grepl("^EXCISED", age))) {
	  warning("No excised fragment for ", filename, ", returning NULL")
	  return(NULL)
  }

  # Extracting the lines pertaining to the reference and assembled contig
  ref    <- grep("^First  seq", age, value = TRUE)
  contig <- grep("^Second seq", age, value = TRUE)
  stopifnot(length(ref) == length(contig))


  # Extracting some info from those lines
  # Sequence name (from fasta) is truncated at blank so it can be parsed by fasta
  ref_start <- as.numeric(sub("^.*\\[\\s*(\\d+),.*", "\\1", ref))
  ref_end   <- as.numeric(sub("^.*\\[\\s*\\d+,\\s*(\\d+)\\].*", "\\1", ref))
  ref_len   <- as.numeric(sub("^.*=>\\s*(\\d+)\\s+nucs.*", "\\1", ref))
  ref_name  <- sub(".*'([[:graph:]]+).*'", "\\1", ref)

  contig_start <- as.numeric(sub("^.*\\[\\s*(\\d+),.*", "\\1", contig))
  contig_end   <- as.numeric(sub("^.*\\[\\s*\\d+,\\s*(\\d+)\\].*", "\\1", contig))
  contig_len   <- as.numeric(sub("^.*=>\\s*(\\d+)\\s+nucs.*", "\\1", contig))
  contig_name  <- sub(".*'([[:graph:]]+).*'", "\\1", contig)

  # Extracting some metadata
  score <- grep("^Score", age, value = TRUE)
  stopifnot(length(score) == 1)
  score <- as.numeric(regmatches(score, regexpr("\\d+", score)))

  identity <- grep("^Identic", age, value = TRUE)
  stopifnot(length(identity) == 1)
  identity <- as.numeric(sub("^.*\\(\\s*(\\d+)%\\)\\s*nucs.*", "\\1", identity))

  gaps <- grep("^Gaps", age, value = TRUE)
  stopifnot(length(gaps) == 1)
  gaps <- as.numeric(sub("^.*\\(\\s*(\\d+)%\\)\\s*nucs.*", "\\1", gaps))

  time <- grep("^Alignment time is", age, value = TRUE)
  stopifnot(length(time) == 1)
  time <- as.numeric(regmatches(time, regexpr("\\d+\\.?\\d*", time)))

  # We have to delimit the lines in which to search for excised regions
  ex_range <- c(grep("^EXCISED REGION\\(S\\):$", age), grep("^Identity at breakpoints", age))
  stopifnot(length(ex_range) == 2 && diff(ex_range) >= 1)
  ex_lines <- age[ex_range[1]:ex_range[2]]

  # Now we initialize a list() that will contain the info on excised regions
  excised <- list()

  # We extract the lines containing the excision info
  ref_excised    <- grep("^ first  seq", ex_lines, value = TRUE)
  contig_excised <- grep("^ second seq", ex_lines, value = TRUE)
  stopifnot(length(ref_excised) >= 1 && length(ref_excised) == length(contig_excised))

  for(i in 1:length(ref_excised)) {
    # Initializing the sub-element of the list
    excised[[i]] <- list()
    excised[[i]]$ref$len <- as.numeric(sub("^.*\\s+(\\d+)\\s+nucs.*", "\\1", ref_excised[i]))
    excised[[i]]$ref$start <- as.numeric(sub("^.*nucs\\s+\\[\\s*(\\d+),.*", "\\1", ref_excised[i]))
    excised[[i]]$ref$end <- as.numeric(sub("^.*\\[.*,(\\d+)\\s*\\].*", "\\1", ref_excised[i]))

    excised[[i]]$contig$len <- as.numeric(sub("^.*\\s+(\\d+)\\s+nucs.*", "\\1", contig_excised[i]))
    excised[[i]]$contig$start <- as.numeric(sub("^.*nucs\\s+\\[\\s*(\\d+),.*", "\\1", contig_excised[i]))
    excised[[i]]$contig$end <- as.numeric(sub("^.*\\[.*,(\\d+)\\s*\\].*", "\\1", contig_excised[i]))
  }

  # Preparing the list object to be returned
  output <- list(reference = list(name = ref_name, range = c(ref_start, ref_end), length = ref_len),
		 contig = list(name = contig_name, range = c(contig_start, contig_end), length = contig_len),
		 score = score,
		 identity = identity,
		 gaps = gaps,
		 time = time,
		 excised = excised)

  return(output)
}

