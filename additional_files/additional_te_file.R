#!/prg/R/4.0/bin/Rscript

# Making a supplementary CSV file from the polymorphic TEs metadata in te_analysis/polymorphic_tes.tsv

# DEPENDENCY: te_analysis/polymorphic_tes.tsv
polymorphic_tes <- read.table("../te_analysis/polymorphic_tes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Reordering the columns as we want them in the output
polymorphic_tes <- polymorphic_tes[, c("qseqid", "chrom", "start", "end", "svtype", "qlen", "alt_freq", "name", "slen", 
				       "length", "evalue", "pident", "te_group", "class", "subclass", "order", "superfamily", 
				       "family", "description")]

names(polymorphic_tes) <- c("sv_id", "sv_chrom", "sv_start_pos", "sv_end_pos", "sv_type", "sv_length", "sv_population_frequency", 
			    "soytedb_blastn_match", "soytedb_te_length", "blastn_alignment_length", "blastn_evalue", 
			    "blastn_percent_identity", "te_group", "te_class", "te_subclass", "te_order", "te_superfamily", 
			    "te_family", "te_description")

# OUTPUT : additional_files/additional_te_file.csv
write.table(polymorphic_tes, file = "additional_te_file.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

