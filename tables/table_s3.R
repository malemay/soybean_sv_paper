#!/prg/R/4.0/bin/Rscript

# Formatting the data for Table S3
# DEPENDENCY : tables/table_s3_data.txt
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
# OUTPUT : tables/table_s3.csv
write.table(table_s3, file = "table_s3.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

