# This file is kept here only for compatibility with the original Makefile
# Metainfo addition is now done within the breakpoint refinement pipeline
# Therefore, this code only copies the files ending with "_realigned.vcf"
# to files ending with "_realigned_metainfo.vcf"

# Adding the metainfo for all samples in turn
# DEPENDENCY : utilities/line_ids.txt
for(i in read.table("../utilities/line_ids.txt", header = TRUE, stringsAsFactors = FALSE)[[1]]) {
# DEPENDENCY : nanopore_sv_calling/AC2001/AC2001_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/ALTA/ALTA_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/MAPLE_ISLE/MAPLE_ISLE_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/MAPLE_PRESTO/MAPLE_PRESTO_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_09_35C/OAC_09_35C_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_CARMAN/OAC_CARMAN_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_DRAYTON/OAC_DRAYTON_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_EMBRO/OAC_EMBRO_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_LAKEVIEW/OAC_LAKEVIEW_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_MADOC/OAC_MADOC_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_OXFORD/OAC_OXFORD_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_PETREL_2/OAC_PETREL_2_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_PRUDENCE/OAC_PRUDENCE_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_STRATFORD/OAC_STRATFORD_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/OT09-03/OT09-03_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/QS5091.50J/QS5091.50J_realigned.vcf
# DEPENDENCY : nanopore_sv_calling/ROLAND/ROLAND_realigned.vcf
	input_file = paste0(i, "/", i, "_realigned.vcf")

	file.copy(input_file, paste0(i, "_realigned_metainfo.vcf"))
}

