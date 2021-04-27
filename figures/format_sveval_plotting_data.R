# This function takes the output of sveval and prepares a data.frame
# for plotting the proportion of true positive and false positive
# calls out of it. It allows selecting either DEL, INS, INV or DUP
format_sveval_plotting_data <- function(x, svtype) {
	# First extracting the relevant SV type
	x <- x[[svtype]]

	# Then extracting the true positive and false positive calls, and only the required columns
	tp <- x$TP
	tp <- tp[, c("size", "ClusterIDs")]
	fp <- x$FP
	fp <- fp[, c("size", "ClusterIDs")]

	# Then we add a column indicating which are true and false positives and merge them
	tp$truepos <- TRUE
	fp$truepos <- FALSE
	gr <- c(tp, fp)

	# The metadata columns are all we need so we coerce them to a data.frame
	plotting_data <- as.data.frame(mcols(gr))

	# Now the only thing left to do is format the ClusterIDs column appropriately
	# We need to:
	# 1- Extract only the name of the program
	# 2- Make sure that any program is mentioned only once
	# 3- Sort the program names alphabetically
	# 4- Paste them back into a single string
	# 5- Coerce them as factor by decreasing number of calls per combination
	for(i in 1:nrow(plotting_data))  {
		i_data <- plotting_data[i, "ClusterIDs"]
		i_data <- strsplit(i_data, ":")[[1]]
		i_data <- sapply(strsplit(i_data, "_"), function(x) x[[1]])
		i_data <- sort(unique(i_data))
		plotting_data[i, "ClusterIDs"] <- paste0(i_data, collapse = ":")
	}

	# Getting the factor levels from the sorting the combinations by number of calls
	lvls <- names(sort(table(plotting_data$ClusterIDs), decreasing = TRUE))
	plotting_data$ClusterIDs <- factor(plotting_data$ClusterIDs, levels = lvls)

	return(plotting_data)
}

