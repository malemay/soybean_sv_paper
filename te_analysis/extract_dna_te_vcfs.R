#!/prg/R/4.0/bin/Rscript

# This code extracts the vcfs of the SVs representing TEs in the database that
# were matched by at least three SVs in our dataset

# Loading the required package
library(VariantAnnotation)

# DEPENDENCY : query_all.vcf
vcf <- readVcf("query_all.vcf", param = ScanVcfParam(info = c("SVTYPE", "AF"), geno = NA))

# DEPENDENCY : polymorphic_tes.tsv
polymorphic_tes <- read.table("polymorphic_tes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extracting a character vector of DNA TE families that are highly polymorphic (their name occurs >= 3 in the polymorphic_tes)
polymorphic_dna_tes <- table(grep("^D", polymorphic_tes$name, value = TRUE))
polymorphic_dna_tes <- polymorphic_dna_tes[polymorphic_dna_tes >= 3]

# Outputting a single vcf file for each of the DNA families in polymorphic_dna_tes
for(i in names(polymorphic_dna_tes)) {
	i_name <- gsub("-", "_", i)
	ids <- polymorphic_tes[polymorphic_tes$name == i, "qseqid"]
	i_vcf <- vcf[rownames(vcf) %in% ids]
	writeVcf(i_vcf, filename = paste0("multiple_alignments/", i_name, ".vcf"))
}

