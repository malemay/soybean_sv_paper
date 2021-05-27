#!/prg/R/4.0/bin/Rscript

# Loading the dataset about breakpoint refinement
# DEPENDENCY : nanopore_sv_calling/all_metainfo.RData
load("../nanopore_sv_calling/all_metainfo.RData")

# Generating a table of the number of (not) realigned variants for each sample
bp_table <- table(all_metainfo$flag, all_metainfo$sample)
bp_table <- t(as.matrix(bp_table))

bp_table <- data.frame(sample = rownames(bp_table),
		       refined = bp_table[, "REALIGNED"],
		       belowthr = bp_table[, "BELOW_TRESHOLDS"],
		       svlen = bp_table[, "MAX_SVLEN"],
		       noalignment = bp_table[, "NO_ALIGNMENT"],
		       stringsAsFactors = FALSE)

# Removing OAC Petrel and formatting the names of the samples
bp_table <- bp_table[bp_table$sample != "OAC_PETREL", ]

line_names <- c("MAPLE_PRESTO" = "Maple Presto",
		"OAC_CARMAN" = "OAC Carman",
		"OAC_EMBRO" = "OAC Embro",
		"QS5091.50J" = "QS5091.50j")

bp_table$sample <- line_names[bp_table$sample]

# Writing to file for inclusion in the .tex document
# OUTPUT : tables/table_s4.csv
write.table(bp_table, file = "table_s4.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

