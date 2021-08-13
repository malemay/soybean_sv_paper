#!/prg/R/4.0/bin/Rscript

# Reading in the csv file with the full table
# DEPENDENCY : tables/lab_methods_table.csv
lab_methods <- read.table("lab_methods_table.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Preparing the second table, which will contain metadata on the sequencing runs
table_s8 <- lab_methods[, c(1:2, 9:13)]
names(table_s8) <- c("sample", "date", "pores", "mass", "runtime", "yield", "nfifty")
table_s8$nfifty <- sprintf("%.1f", table_s8$nfifty)

# Writing to file
# OUTPUT : tables/table_s8.csv
write.table(table_s8, file = "table_s8.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

