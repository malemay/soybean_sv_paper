#!/usr/bin/Rscript

# Loading the datasets
# DEPENDENCY : raw_variants/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData
load("/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_nanopore/tests/raw_variants/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData")
raw_variants_rates <- sveval_nogeno_rates
names(raw_variants_rates) <- c("DEL", "INS")

# DEPENDENCY : realigned_variants/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData
load("/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_nanopore/tests/realigned_variants/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData")
realigned_variants_rates <- sveval_nogeno_rates
names(realigned_variants_rates) <- c("DEL", "INS")

rm(sveval_nogeno_rates)

# Sourcing the function used to prepare the data for plotting
# DEPENDENCY : make_plot_data.R
source("~/scripts/make_plot_data.R")

# Plot-ready data.frames for deletions and insertions are prepared
raw_deletions <- make_plot_data(raw_variants_rates, "DEL")
raw_deletions$type <- "raw"
realigned_deletions <- make_plot_data(realigned_variants_rates, "DEL")
realigned_deletions$type <- "realigned"
deletions <- rbind(raw_deletions, realigned_deletions)
deletions$linegroup2 <- paste0(deletions$linegroup, deletions$type)

raw_insertions <- make_plot_data(raw_variants_rates, "INS")
raw_insertions$type <- "raw"
realigned_insertions <- make_plot_data(realigned_variants_rates, "INS")
realigned_insertions$type <- "realigned"
insertions <- rbind(raw_insertions, realigned_insertions)
insertions$linegroup2 <- paste0(insertions$linegroup, insertions$type)

# Removing OAC PETREL from the data and renaming the cultivars
deletions <- deletions[deletions$cultivar != "CAD1022", ]
insertions <- insertions[insertions$cultivar != "CAD1022", ]

line_names <- c("CAD1010" = "Maple Presto",
		"CAD1052" = "OAC Embro",
		"CAD1064" = "OAC Carman",
		"CAD1070" = "QS5091.50j")

deletions$cultivar <- line_names[deletions$cultivar]
insertions$cultivar <- line_names[insertions$cultivar]

# Removing the "[30-50[" and "all" size classes
deletions <- deletions[!deletions$size_class %in% c("[30-50[", "all"), ]
deletions$size_class <- droplevels(deletions$size_class)

insertions <- insertions[!insertions$size_class %in% c("[30-50[", "all"), ]
insertions$size_class <- droplevels(insertions$size_class)

# Saving to file
save(deletions, file = "deletions.RData")
save(insertions, file = "insertions.RData")

