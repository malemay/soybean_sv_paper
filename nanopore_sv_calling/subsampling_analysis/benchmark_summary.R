#!/prg/R/4.0/bin/Rscript

# Looping over all the RDS files to retrieve the results
rates_list <- list()

# DEPENDENCY : nanopore_sv_calling/subsampling_analysis/SUBSAMPLE_BENCHMARKS
rds_files <- dir("benchmark_files/", pattern = "\\.rds$")

for(i in rds_files) {
	rates_list[[i]] <- readRDS(paste0("benchmark_files/", i))
}

# Merging the outputs of all analyses together
deletions  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "DEL")))
insertions  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "INS")))
inversions  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "INV")))
duplications  <- do.call("rbind", lapply(rates_list, function(x) `[[`(x, "DUP")))

# Adding some columns and merging the datasets
deletions$svtype <- "DEL"
insertions$svtype <- "INS"
inversions$svtype <- "INV"
duplications$svtype <- "DUP"

benchmark_rates <- rbind(deletions, insertions, inversions, duplications)

# Filtering out the 30-50 size class
benchmark_rates <- benchmark_rates[benchmark_rates$size_class != "[30-50[", ]

# Checking if all threshold values yield the same result
thresholds_identical <- tapply(X = benchmark_rates$sensitivity, 
			       INDEX = list(benchmark_rates$cultivar, benchmark_rates$size_class, benchmark_rates$svtype),
			       FUN = function(x) all(x == x[1]))
stopifnot(all(thresholds_identical))

# They are, so we can keep only the value where the threshold equals 0
benchmark_rates <- benchmark_rates[benchmark_rates$threshold == 0, ]

# Also removing some un-neded columns and row names
benchmark_rates$threshold <- NULL
benchmark_rates$pipeline <- NULL
benchmark_rates$precision_shrunk <- NULL
rownames(benchmark_rates) <- NULL

# Adding some columns and reformatting the "cultivars" column
benchmark_rates$frac <- as.numeric(sub(".*frac(.*)_seed.*", "\\1", benchmark_rates$cultivar))
benchmark_rates$rep <- as.numeric(sub(".*seed(.*)", "\\1", benchmark_rates$cultivar))
benchmark_rates$cultivar <- sub("(.*)_frac.*", "\\1", benchmark_rates$cultivar)
benchmark_rates$rep <- 
	unlist(tapply(benchmark_rates$rep, list(benchmark_rates$cultivar, benchmark_rates$frac), function(x) as.integer(as.factor(x))))

# Loading the sequencing depth per sample so we can translate the fraction into actual sequencing depths
# DEPENDENCY : depth_distributions/average_depth.RData
load("../../depth_distributions/average_depth.RData")
depth <- average_depth$depth
names(depth) <- average_depth$sample
benchmark_rates$depth <- benchmark_rates$frac * depth[benchmark_rates$cultivar]

# Displaying the cases with NA values
benchmark_rates[!complete.cases(benchmark_rates),] 

# Formatting the benchmark data for publication-ready plotting
# Each sample/sequencing depth combination will be plotted using
#  a dot for the median and errorbars for the minimum and maximum 
#  of 5 replicates
benchmark_plotting <- unique(benchmark_rates[, c("size_class", "cultivar", "svtype", "frac", "depth")])

# Creating columns to be filled
benchmark_plotting$s_median <- numeric(nrow(benchmark_plotting))
benchmark_plotting$s_min <- numeric(nrow(benchmark_plotting))
benchmark_plotting$s_max <- numeric(nrow(benchmark_plotting))
benchmark_plotting$p_median <- numeric(nrow(benchmark_plotting))
benchmark_plotting$p_min <- numeric(nrow(benchmark_plotting))
benchmark_plotting$p_max <- numeric(nrow(benchmark_plotting))

# Looping over the rows to fill the values
for(i in 1:nrow(benchmark_plotting)) {
	tmp_df <- benchmark_rates[benchmark_rates$size_class == benchmark_plotting[i, "size_class"] &
				  benchmark_rates$cultivar == benchmark_plotting[i, "cultivar"] &
				  benchmark_rates$svtype == benchmark_plotting[i, "svtype"] &
				  benchmark_rates$frac == benchmark_plotting[i, "frac"], ]

	stopifnot(nrow(tmp_df) == 5)

	benchmark_plotting[i, "s_median"] <- median(tmp_df$sensitivity, na.rm = TRUE)
	benchmark_plotting[i, "s_min"] <- min(tmp_df$sensitivity, na.rm = TRUE)
	benchmark_plotting[i, "s_max"] <- max(tmp_df$sensitivity, na.rm = TRUE)
	benchmark_plotting[i, "p_median"] <- median(tmp_df$precision, na.rm = TRUE)
	benchmark_plotting[i, "p_min"] <- min(tmp_df$precision, na.rm = TRUE)
	benchmark_plotting[i, "p_max"] <- max(tmp_df$precision, na.rm = TRUE)
}

# We need to add sensitivity and precision values of 1 for each sample at its full sequencing depth
full_rates <- unique(benchmark_plotting[, c("size_class", "cultivar", "svtype")])
full_rates$frac <- 1
full_rates$depth <- depth[full_rates$cultivar]
full_rates$s_median  <- full_rates$s_min <- full_rates$s_max <- full_rates$p_median <- full_rates$p_min <- full_rates$p_max <- 1

# We add those to the benchmark_plotting data.frame
benchmark_plotting <- rbind(benchmark_plotting, full_rates)

# Saving the data to a file for plotting from other scripts
save(benchmark_plotting, file = "benchmark_data.RData")

