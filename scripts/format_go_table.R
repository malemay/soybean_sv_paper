# Writing a function that will be used to format GO/PFAM enrichment tables (Tables S6, S7 and S8)
format_go_table <- function(x){
	x <- x[x$Pvalue < 0.05, ]
	x$Pvalue <- sprintf("%.4g", x$Pvalue)
	x$OddsRatio <- sprintf("%.2f", x$OddsRatio)
	x$ExpCount <- sprintf("%.2f", x$ExpCount)
	x
}

