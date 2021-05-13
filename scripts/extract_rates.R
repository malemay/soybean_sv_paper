# Writing a function which generates a data.frame of precision
#  and sensitivity rates. This time, working directly from the output containing
#  the results of the evaluation at all quality thresholds

# IMPORTANT : The version of sveval pacakge currently installed has been modified
# to allow for the output of all quality thresholds to be returned. This in turns
# allows this function to extract the precision/sensitivity rates for all
# quality thresholds and various size classes.
extract_rates <- function(x, svtypes = c("DEL", "INS"), cultivar = NULL, pipeline = NULL) {

	# Checking that the sample name is provided
	if(is.null(cultivar) || is.null(pipeline)) stop("Cultivar name and pipeline must both be provided")

	# Creating a template data.frame to be filled as the loop goes forward
	template <- data.frame(size_class = character(),
			       sensitivity = numeric(),
			       precision = numeric(),
			       precision_shrunk = numeric(),
			       pipeline = character(),
			       cultivar = character(),
			       threshold = numeric(),
			       stringsAsFactors = FALSE)

	# Creating a dataframe of labels and size limits to iterate on
	sizes <- data.frame(label = c("all", "[30-50[", "[50-100[", "[100-1000[", "[1000-10000[", "[10000+["),
			    lower = c(0, 30, 50, 100, 1000, 10000),
			    upper = c(Inf, 50, 100, 1000, 10000, Inf),
			    stringsAsFactors = FALSE)

	# Initializing a list that will contain the output
	# It is a list with data.frames for each type of SV
	output <- list()

	# Looping over the sv types to be studied
	for(sv_type in svtypes) {

		# Initializing a data.frame for that SV type
		output[[sv_type]] <- template


		# Creating a loop to iterate over the sizes
		for(i in 1:nrow(sizes)) {

			# Extracting the information from the data.frame for this iteration
			label <- sizes[i, "label"]
			lower <- sizes[i, "lower"]
			upper <- sizes[i, "upper"]

			# Getting the number of true positive events
			tp <- sapply(x$eval, function(x) sum(x$regions[[sv_type]]$TP$size >= lower & x$regions[[sv_type]]$TP$size < upper))
			# Getting the number of true positive baseline events
			tpb <- sapply(x$eval, function(x) sum(x$regions[[sv_type]]$TP.baseline$size >= lower & x$regions[[sv_type]]$TP.baseline$size < upper))
			# Getting the number of false positive events
			fp <- sapply(x$eval, function(x) sum(x$regions[[sv_type]]$FP$size >= lower & x$regions[[sv_type]]$FP$size < upper))
			# Getting the number of false negative events
			fn <- sapply(x$eval, function(x) sum(x$regions[[sv_type]]$FN$size >= lower & x$regions[[sv_type]]$FN$size < upper))

			# Creating the deletion data.frame row for this iteration
			i_df <-  data.frame(size_class = label,
					    sensitivity = tpb / (tpb + fn),
					    precision = tp / (tp + fp),
					    precision_shrunk = tpb / (tpb + fp),
					    pipeline = pipeline,
					    cultivar = cultivar,
					    threshold = x$qual_ths,
					    stringsAsFactors = FALSE)

			# Adding the data.frame rows thus generated to the main data.frames
			output[[sv_type]]  <- rbind(output[[sv_type]], i_df)
		}
	}

	return(output)
}

