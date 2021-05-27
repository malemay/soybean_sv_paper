#!/prg/R/4.0/bin/Rscript

# Reading in the csv file with the full table
lab_methods <- read.table("lab_methods_table.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Preparing the second table, which will contain metadata on the sequencing runs
table_s2 <- lab_methods[, c(1:2, 9:13)]
names(table_s2) <- c("sample", "date", "pores", "mass", "runtime", "yield", "nfifty")
table_s2$nfifty <- sprintf("%.1f", table_s2$nfifty)

# Writing to file
write.table(table_s2, file = "table_s2.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

