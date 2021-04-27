#!/usr/bin/Rscript

# In this script, I prepare three circos data input files that 
# will be used to highlight regions of interest in the lineplot
# and histogram tracks

# For the SV histogram track, I will highlight the regions with
# the 10% highest numbers of SVs (deletions and insertions
# combined)

# For the TE lineplot tracks, I will highlight the regions
# with the 10% highest ratios of number of polymoprhic TEs
# to the number of reference TEs

# Reading the SV data track
sv_track <- read.table("sv_counts.txt", header = FALSE, stringsAsFactors = FALSE)
sv_track[[5]] <- apply(sv_track, 1, function(x) sum(as.numeric(strsplit(x[4], ",")[[1]])))
sv_track[[6]] <- sv_track[[5]] >= quantile(sv_track[[5]], 0.90)

# Writing the bins corresponding to the condition to file
write.table(sv_track[sv_track[[6]], 1:3], col.names = FALSE, row.names = FALSE,
	    quote = FALSE, sep = "\t", file = "sv_highlights.txt")

# Reading the LTR elements tracks
ref_ltr  <- read.table("ref_ltr.txt",  header = FALSE, stringsAsFactors = FALSE)
poly_ltr <- read.table("poly_ltr.txt", header = FALSE, stringsAsFactors = FALSE)

ltr_track <- ref_ltr
ltr_track[[4]] <- poly_ltr[[4]] / ref_ltr[[4]]

ltr_track[[5]] <- ltr_track[[4]] >= quantile(ltr_track[[4]], 0.90)

# Writing the track to file for plotting with circos
write.table(ltr_track[ltr_track[[5]], 1:3], col.names = FALSE, row.names = FALSE,
	    quote = FALSE, sep = "\t", file = "ltr_highlights.txt")



# Reading the DNA elements tracks
ref_dna  <- read.table("ref_dna.txt",  header = FALSE, stringsAsFactors = FALSE)
poly_dna <- read.table("poly_dna.txt", header = FALSE, stringsAsFactors = FALSE)

dna_track <- ref_dna
dna_track[[4]] <- poly_dna[[4]] / ref_dna[[4]]

dna_track[[5]] <- dna_track[[4]] >= quantile(dna_track[[4]], 0.90)

# Writing the track to file for plotting with circos
write.table(dna_track[dna_track[[5]], 1:3], col.names = FALSE, row.names = FALSE,
	    quote = FALSE, sep = "\t", file = "dna_highlights.txt")

