#!/prg/R/4.0/bin/Rscript

# This code formats the supplementary tables S1 and S2 pertaining to lab methods

# Reading in the csv file with the full table
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
write.table(table_s1, file = "table_s1.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")


# Now preparing the second table, which will contain metadata on the sequencing runs
table_s2 <- lab_methods[, c(1:2, 9:13)]
names(table_s2) <- c("sample", "date", "pores", "mass", "runtime", "yield", "nfifty")
table_s2$nfifty <- sprintf("%.1f", table_s2$nfifty)

# Writing to file
write.table(table_s2, file = "table_s2.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")



# Formatting the data for Table S3
table_s3 <- read.table("table_s3_data.txt", header = TRUE, stringsAsFactors = FALSE)

# Reformatting the names of the lines in the sample column
ind_names <- c("AC2001" = "AC2001",
	       "ALTA" = "Alta",
	       "MAPLE_ISLE" = "Maple Isle",
	       "MAPLE_PRESTO" = "Maple Presto",
	       "OAC_09_35C" = "OAC 09-35C",
	       "OAC_CARMAN" = "OAC Carman",
	       "OAC_DRAYTON" = "OAC Drayton",
	       "OAC_EMBRO" = "OAC Embro",
	       "OAC_LAKEVIEW" = "OAC Lakeview",
	       "OAC_MADOC" = "OAC Madoc",
	       "OAC_OXFORD" = "OAC Oxford",
	       "OAC_PETREL_2" = "OAC Petrel",
	       "OAC_PRUDENCE" = "OAC Prudence",
	       "OAC_STRATFORD" = "OAC Stratford",
	       "OT09-03" = "OT09-03",
	       "QS5091.50J" = "QS5091.50j",
	       "ROLAND" = "Roland")

table_s3$sample <- ind_names[table_s3$sample]

# Write to file as a csv
write.table(table_s3, file = "table_s3.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

