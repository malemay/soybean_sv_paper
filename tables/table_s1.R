#!/prg/R/4.0/bin/Rscript

# Reading in the csv file with the full table
# DEPENDENCY : tables/lab_methods_table.csv
lab_methods <- read.table("lab_methods_table.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# The first table will contain the first 8 columns, except date
table_s1 <- lab_methods[, 1:8]
table_s1 <- table_s1[, -2]

# Creating some lookup tables to recode the columns
protocol  <- c("CTAB" = "CTAB", 
	      "Qiagen DNeasy Plant Mini Kit" = "DNeasy",
	      "Qiagen Gentra Puregene Tissue Kit" = "Gentra")

grinding <- c("Liquid nitrogen + Qiagen TissueLyser" = "LN + TL",
	      "Lyophilisation + Qiagen TissueLyser" = "CD + TL",
	      "Liquid nitrogen + mortar and pestle" = "LN + MP")

size_selection <- c("None" = "None",
		    "BluePippin High-Pass Plus > 15kb" = "BP 15kb",
		    "BluePippin > 6kb" = "BP 6kb",
		    "Circulomics SRE" = "SRE")

names(table_s1) <- c("sample", "growing", "grinding", "extraction", "rnase", "size", "combined")

# Changing the contents of the columns of Table S1
table_s1$grinding <- grinding[table_s1$grinding]
table_s1$extraction <- protocol[table_s1$extraction]
table_s1$size <- size_selection[table_s1$size]

# Writing the contents of the table to table_s1.csv
# OUTPUT : tables/table_s1.csv
write.table(table_s1, file = "table_s1.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

