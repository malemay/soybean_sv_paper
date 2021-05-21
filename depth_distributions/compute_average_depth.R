# Computing a data.frame with the average Oxford Nanopore sequencing depth of all 17 samples of interest

# This chunk reads in the data on the genome-wide sequencing depth for all samples
depth_data <- list()

# Iterating over all the samples
# DEPENDENCY : utilities/line_ids.txt
for(i in read.table("../utilities/line_ids.txt", header = TRUE, stringsAsFactors = FALSE)[[1]]) {
	# Reading the data as a data.frame
	depth_data[[i]] <- read.table(paste0(i, "_depth.txt"))
	# Ordering the data.frame by number of reads covering the position
	depth_data[[i]] <- depth_data[[i]][order(depth_data[[i]][[1]]), ]
	# Adding a column with the cumulative sum of number of positions covered
	depth_data[[i]][[3]] <- cumsum(depth_data[[i]][[2]])
	depth_data[[i]][[4]] <- depth_data[[i]][[3]] / sum(depth_data[[i]][[2]])
	depth_data[[i]][[5]] <- i
}

# Binding them all together in a single data.frame
depth_data <- do.call("rbind", depth_data)
names(depth_data) <- c("depth", "n", "cumsum", "cumsum_prop", "sample")

# Computing the average sequencing depth for each sample
average_depth <- numeric(length(unique(depth_data$sample)))
names(average_depth) <- unique(depth_data$sample)

# Iterating over the samples
for(i in names(average_depth)) {
	i_data <- depth_data[depth_data$sample == i, ]
	average_depth[i] <- sum(i_data$depth * i_data$n) / sum(i_data$n)
}

# Turning this into a data.frame
average_depth <- data.frame(sample = names(average_depth),
			    depth = unname(average_depth),
			    stringsAsFactors = FALSE)

# Saving the data.frame to file
# OUTPUT : depth_distributions/average_depth.RData
save(average_depth, file = "average_depth.RData")

