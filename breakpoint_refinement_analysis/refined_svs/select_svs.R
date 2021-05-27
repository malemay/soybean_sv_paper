#!/prg/R/4.0/bin/Rscript

# Using a custom R function to re-select the SVs clustered by SVmerge by preferentially
# choosing realigned variants
# DEPENDENCY : scripts/merge_realigned_variants.R
source("../../scripts/merge_realigned_variants.R")

# Lauching the command
merge_aligned_variants("svmerged_preliminary.clustered.vcf",
		       paste0("../../nanopore_sv_calling/", c("OAC_CARMAN", "OAC_EMBRO", "MAPLE_PRESTO", "QS5091.50J"), "_normalized_ids.vcf"),
		       "svmerged.clustered.vcf")

