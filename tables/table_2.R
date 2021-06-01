#!/prg/R/4.0/bin/Rscript

# Loading the rtracklayer package for importing the GFF3 file
library(rtracklayer)

# Printing a table of the number of matches to different TE types according to SV type
# First loading the TE dataset

# DEPENDENCY : te_analysis/polymorphic_tes.tsv
tes <- read.table("../te_analysis/polymorphic_tes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

te_table <- as.matrix(table(tes$te_group, tes$svtype))

te_table <- data.frame(DEL = te_table[, "DEL"],
		       INS = te_table[, "INS"])

# Adding percentage columns to the data.frame
te_table[, "DELprop"] <- te_table[, "DEL"] / sum(te_table[, "DEL"]) * 100
te_table[, "INSprop"] <- te_table[, "INS"] / sum(te_table[, "INS"]) * 100
# Rounding to the first decimal
te_table[, "DELprop"] <- round(te_table[, "DELprop"], 1)
te_table[, "INSprop"] <- round(te_table[, "INSprop"], 1)

# Adding columns with the number of base pairs
te_ins <- tes[tes$svtype == "INS", ]
te_del <- tes[tes$svtype == "DEL", ]

te_table[, "DELkb"] <- round(as.numeric(tapply(te_del$qlen, te_del$te_group, sum)) / 1000)
te_table[, "INSkb"] <- round(as.numeric(tapply(te_ins$qlen, te_ins$te_group, sum)) / 1000)

te_table[, "DELkbprop"] <- te_table[, "DELkb"] / sum(te_table[, "DELkb"]) * 100
te_table[, "INSkbprop"] <- te_table[, "INSkb"] / sum(te_table[, "INSkb"]) * 100

te_table[, "DELkbprop"] <- round(te_table[, "DELkbprop"], 1)
te_table[, "INSkbprop"] <- round(te_table[, "INSkbprop"], 1)

# Formatting the TE type column and reordering the rows
te_table$type <- rownames(te_table)
te_types <- c("Copia" = "Copia LTR retrotransposons", "Gypsy" = "Gypsy LTR retrotransposons", 
	      "NON-LTR" = "Non-LTR retrotransposons", "DNA" = "DNA TE")
te_table$type <- te_types[te_table$type]
te_table <- te_table[c(1, 3, 4, 2),]

# We also want to add the information on these transposable elements in the
# reference genome (their absolute number and bases occupied) so we can compare
# them to the values obtained for the polymorphic elements
# For this we need to read the reference TEs from file
# DEPENDENCY : refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3
ref_tes <- import("../refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3")

# Let us have a look at the number of elements in each class
table(ref_tes$class)

# Now let us create a GRanges object for each of the TE classes of interest
copia  <- ref_tes[ref_tes$class %in% c("LTR/Copia")]
gypsy  <- ref_tes[ref_tes$class %in% c("LTR/Gypsy", "LTR/Gypsy?")]
nonltr <- ref_tes[ref_tes$class %in% c("LINE/L1", "LINE/RTE-BovB")]
dna    <- ref_tes[ref_tes$class %in% c("DNA/CACTA", "DNA/CMC-EnSpm", "DNA/hAT-Ac", "DNA/hAT-Tag1", "DNA/hAT-Tip100",
				       "DNA/Helitron", "DNA/MULE", "DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/TcMar-Stowaway")]

# We will restrict oursevles to those >= 100 bp because this was our size limit on polymoprhic TEs
copia <- copia[width(copia) >= 100]
gypsy <- gypsy[width(gypsy) >= 100]
nonltr <- nonltr[width(nonltr) >= 100]
dna <- dna[width(dna) >= 100]

# Using these GRanges to add the relevant columns to our data.frame
te_table$ref   <- 0
te_table$refkb <- 0

# Adding the number of elements of the given types in the reference
te_table["Copia", "ref"] <- length(copia)
te_table["Gypsy", "ref"] <- length(gypsy)
te_table["NON-LTR", "ref"] <- length(nonltr)
te_table["DNA", "ref"] <- length(dna)

# Adding the total length (in kb) of elements of the given types in the reference
te_table["Copia", "refkb"] <- sum(width(copia)) / 1000
te_table["Gypsy", "refkb"] <- sum(width(gypsy)) / 1000
te_table["NON-LTR", "refkb"] <- sum(width(nonltr)) / 1000
te_table["DNA", "refkb"] <- sum(width(dna)) / 1000
# Adding a column with the total lengths in Mb because the numbers are so large
te_table$refmb <- round(te_table$refkb / 1000)

# Computing columns with the proportions
te_table$refprop <- sprintf("%.1f", te_table$ref / sum(te_table$ref) * 100)
te_table$refmbprop <- sprintf("%.1f", te_table$refmb / sum(te_table$refmb) * 100)

# Writing the table to a .csv file so it can be used by Latex
# OUTPUT : tables/table_2.csv
write.table(te_table, file = "table_2.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

