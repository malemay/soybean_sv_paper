#!/usr/bin/Rscript

# This file prepares the axis labels data tracks for adding to the Circos plot
# These values cannot be based solely on min/max values because it would not yield
# "pretty" values. Therefore, I will look at the values individually and choose
# minimum (most likely 0) and maximum values. These values will also need to be set
# within the <plot> blocks of the plots.conf file

# Reading the dataset of reference LTR elements
ref_ltr <- read.table("ref_ltr.txt", header = FALSE, stringsAsFactors = FALSE)
range(ref_ltr[[4]])
# [1]   15 1547
# We can set min and max for this track to 0 and 1600

# Reading the dataset of polymorphic LTR elements
poly_ltr <- read.table("poly_ltr.txt", header = FALSE, stringsAsFactors = FALSE)
range(poly_ltr[[4]])
# [1]  0 57
# We can set min and max for this track to 0 and 60

# Reading the dataset of reference DNA transposable elements
ref_dna <- read.table("ref_dna.txt", header = FALSE, stringsAsFactors = FALSE)
range(ref_dna[[4]])
# [1]  12 612
# We can set min and max for this track to 0 and 650

# Reading the dataset of reference DNA transposable elements
poly_dna <- read.table("poly_dna.txt", header = FALSE, stringsAsFactors = FALSE)
range(poly_dna[[4]])
# [1] 0 9
# We can set the min and max for this track to 0 and 10


# The ranges will be output to 4 different files:
# dna_axis_min.txt : the minimum values for the DNA TE axes
# dna_axis_max.txt : the maximum values for the DNA TE axes
# ltr_axis_min.txt : the minimum values for the LTR TE axes
# ltr_axis_max.txt : the maximum values for the LTR TE axes

# Let us write a function that takes the two values for each file and saves them to disk
write_axis_file <- function(ref_value, poly_value, output_file, dummy_length = 60000000) {
	ref_line <- c("dummy", "0", "0", ref_value, "color=blue")
	ref_line <- paste0(ref_line, collapse = "\t")
	cat(ref_line, "\n", sep = "", file = output_file, append = FALSE)

	poly_line <- c("dummy", dummy_length, dummy_length, poly_value, "color=red")
	poly_line <- paste0(poly_line, collapse = "\t")
	cat(poly_line, "\n", sep = "", file = output_file, append = TRUE)

	return(invisible(NULL))
}

# Changing the scipen option
options(scipen = 12)

# Writing each of the files
write_axis_file(ref_value = 0, poly_value = 0, output_file = "ltr_axis_min.txt")
write_axis_file(ref_value = 1600, poly_value = 60, output_file = "ltr_axis_max.txt")
write_axis_file(ref_value = 0, poly_value = 0, output_file = "dna_axis_min.txt")
write_axis_file(ref_value = 650, poly_value = 10, output_file = "dna_axis_max.txt")

