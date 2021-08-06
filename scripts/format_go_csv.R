# Writing a function that will be used to format GO/PFAM enrichment tables for output to csv
format_go_csv <- function(x, output_file, go_test){
	names(x) <- c("go_id", "p_value", "odds_ratio", "expected_count", "observed_count", "number_of_genes", "description")
	x <- x[, c("go_id", "description", "number_of_genes", "expected_count", "observed_count", "odds_ratio", "p_value")]
	if(!go_test) names(x)[1] <- "pfam_id"
	write.table(x, file = output_file, sep = ",", col.names = TRUE, row.names = FALSE, quote = TRUE)
	invisible(NULL)
}

