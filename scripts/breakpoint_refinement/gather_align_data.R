# File initially created on Tuesday, September 29, 2020

# This function takes the output of a call to parse_svinfo and
# that of a call to parse_age for the same SV and formats
# the data to be used for comparing the two ranges of coordinates
# into a data.frame

gather_align_data <- function(svinfo, age, window_range, contig_fasta) {

  # We need to extract the data from age in case more than one alignment was found
  excise_ref_len <- sapply(age$excise, function(x) x$ref$len)
  excise_ref_start <- sapply(age$excise, function(x) x$ref$start)
  excise_ref_end <- sapply(age$excise, function(x) x$ref$end)

  excise_contig_len <- sapply(age$excise, function(x) x$contig$len)
  excise_contig_start <- sapply(age$excise, function(x) x$contig$start)
  excise_contig_end <- sapply(age$excise, function(x) x$contig$end)
  
  # Checking that the input is suitable
  if(length(excise_ref_len)      < 1) return(NULL)
  if(length(excise_ref_start)    < 1) return(NULL)
  if(length(excise_ref_end)      < 1) return(NULL)
  if(length(excise_contig_len)   < 1) return(NULL)
  if(length(excise_contig_start) < 1) return(NULL)
  if(length(excise_contig_end)   < 1) return(NULL)

  output <- data.frame(chr = svinfo$chr,
		       svtype = svinfo$svtype,
  		       start = svinfo$start,
		       svlen = svinfo$svlen,
		       alt = svinfo$alt,
		       imprecise = svinfo$imprecise,
		       age_refname = age$reference$name,
		       age_refstart = age$reference$range[1],
		       age_refend = age$reference$range[2],
		       age_reflen = age$reference$length,
		       age_contigname = age$contig$name,
		       age_contigstart = age$contig$range[1],
		       age_contigend = age$contig$range[2],
		       age_contiglen = age$contig$length,
		       score = age$score,
		       identity = age$identity,
		       gaps = age$gaps,
		       time = age$time,
		       excise_ref_len = excise_ref_len,
		       excise_ref_start = excise_ref_start,
		       excise_ref_end = excise_ref_end,
		       excise_contig_len = excise_contig_len,
		       excise_contig_start = excise_contig_start,
		       excise_contig_end = excise_contig_end,
		       stringsAsFactors = FALSE)

  # Adding a series of variables to make the SVs truly comparable
  output$age_start  <- output$excise_ref_start - window_range + output$start
  output$age_length <- abs(excise_contig_len - excise_ref_len)

  # Getting the alternate sequence of insertions from the contig fasta
  if(all(output$svtype == "INS")) {
    ins_range <- GenomicRanges::GRanges(seqnames = output$age_contigname,
					IRanges::IRanges(start = pmin(output$excise_contig_start, 
							     output$excise_contig_end),
						end = pmax(output$excise_contig_start,
							   output$excise_contig_end)))
    Rsamtools::indexFa(contig_fasta)
    output$age_alt <- as.character(Rsamtools::scanFa(contig_fasta, ins_range))

    # We have to set it to reverse complement if it was aligned in reverse
    for(i in 1:nrow(output)) {
      if(output[i, "age_contigstart"] > output[i, "age_contigend"]) {
	output[i, "age_alt"] <- revcomp(output[i, "age_alt"])
      }

    # We also compute the Levenshtein distance
    output[i, "lev_dist"] <- adist(output[i, "alt"], output[i, "age_alt"])
    }

  # Then we compute the distance ratio as the Levensthein distance divided by the largest of the two sequences
    output$dist_ratio <- output$lev_dist / pmax(output$svlen, output$excise_contig_len)

    # We also determine the distance between the two insertion points
    output$offset <- abs(output$start - output$age_start)

  } else {
    output$age_alt <- NA
    output$lev_dist <- NA
    output$dist_ratio <- NA
    output$offset <- NA
  }

  # For deletions we want to compute the reciprocal overlap between the two deletions
  if(all(output$svtype == "DEL")) {
	  # Creating IRanges objects for each of the two sets
	  sniffles_ranges <- IRanges::IRanges(start = output$start, width = output$svlen)
	  age_ranges      <- IRanges::IRanges(start = output$age_start, width = output$excise_ref_len)

	  # Computing the overlapping ranges between the two
	  ol_ranges <- IRanges::pintersect(sniffles_ranges, age_ranges, resolve.empty = "start.x")
	  # Computing a column for the reciprocal overlap (minimum of the two overlaps)
	  output$r_overlap <- pmin(IRanges::width(ol_ranges) / IRanges::width(sniffles_ranges), 
				   IRanges::width(ol_ranges) / IRanges::width(age_ranges))
	  
	  # If there is no deletion, we set the overlap to 0 (otherwise it is NaN)
	  output$r_overlap[IRanges::width(age_ranges) == 0] <- 0

  } else {
	  output$r_overlap <- NA
  }
  
  return(output)
}

