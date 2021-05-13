# File initially created on September 11, 2020

# This function takes a list of two data.frames (one for deletions, the other for insertions)
#  containing structural variant comparison benchmarks and processes the data in order
#  to make it suitable for plotting with ggplot2

make_plot_data <- function(x, type) {

	# Making adjustments for backwards compatibility
	if(type == "insertions") type <- "INS"
	if(type == "deletions")  type <- "DEL"

	# Checking inputs
	stopifnot(type %in% c("DEL", "INS", "INV", "DUP"))

	# Extracting the deletion or insertion data depending on input
	x <- x[[type]]

	# Creating a variable to group the lines together
	x$linegroup <- paste0(x$pipeline, x$cultivar)

	# Ordering the levels of size_class
	x$size_class <- factor(x$size_class, levels = c("[30-50[", "[50-100[", "[100-1000[", 
							"[1000-10000[", "[10000+[", "all"))

	# Returning the results, ready for plotting
	return(x)
}

