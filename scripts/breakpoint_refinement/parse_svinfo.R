# File initially created on Tuesday, September 29, 2020

# A function that parses a line in a VCF file in order to extract metadata
# about the structural variants it represents

parse_svinfo <- function(vcf_line) {

  stopifnot(length(vcf_line) == 1)

  # Splitting the vcf line into its fields
  split_line <- strsplit(vcf_line, "\\s+")[[1]]

  # Assigning some elements directly
  chr <- split_line[1]
  start <- as.numeric(split_line[2])
  alt <- split_line[5]

  # Some require a bit of processing
  svtype <- sub(".*SVTYPE=([[:upper:]]+).*", "\\1", split_line[8])
  svlen  <- as.numeric(sub(".*SVLEN=-?(\\d+).*", "\\1", split_line[8]))
  imprecise <- grepl("IMPRECISE", split_line[8])

  list(chr = chr, svtype = svtype, start = start, svlen = svlen, alt = alt, imprecise = imprecise)
}

