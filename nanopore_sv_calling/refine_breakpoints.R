# Loading the necessary functions
# DEPENDENCY : scripts/breakpoint_refinement/breakpoint_refinement.R
source("../scripts/breakpoint_refinement/breakpoint_refinement.R")

# Reading the location of the executables from the command line
executables <- commandArgs(trailingOnly = TRUE)
samtools <- executables[1]
minimap2 <- executables[2]
age <- executables[3]
wtdbg2 <- executables[4]
wtpoa_cns <- executables[5]

# DEPENDENCY : utilities/line_ids.txt
for (i in read.table("../utilities/line_ids.txt", header = TRUE, stringsAsFactors = FALSE)[[1]]) {
	dir.create(i)
	setwd(i)
# DEPENDENCY : nanopore_sv_calling/AC2001_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/ALTA_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/MAPLE_ISLE_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/MAPLE_PRESTO_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_09_35C_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_CARMAN_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_DRAYTON_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_EMBRO_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_LAKEVIEW_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_MADOC_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_OXFORD_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_PETREL_2_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_PRUDENCE_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OAC_STRATFORD_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/OT09-03_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/QS5091.50J_hom70_filtered.vcf
# DEPENDENCY : nanopore_sv_calling/ROLAND_hom70_filtered.vcf
	# Refining the breakpoints of the filtered Sniffles vcf using the following parameters
	refine_breakpoints(input_vcf = paste0("../", i, "_hom70_filtered.vcf"),
# OUTPUT : nanopore_sv_calling/AC2001/AC2001_realigned.vcf
# OUTPUT : nanopore_sv_calling/ALTA/ALTA_realigned.vcf
# OUTPUT : nanopore_sv_calling/MAPLE_ISLE/MAPLE_ISLE_realigned.vcf
# OUTPUT : nanopore_sv_calling/MAPLE_PRESTO/MAPLE_PRESTO_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_09_35C/OAC_09_35C_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_CARMAN/OAC_CARMAN_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_DRAYTON/OAC_DRAYTON_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_EMBRO/OAC_EMBRO_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_LAKEVIEW/OAC_LAKEVIEW_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_MADOC/OAC_MADOC_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_OXFORD/OAC_OXFORD_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_PETREL_2/OAC_PETREL_2_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_PRUDENCE/OAC_PRUDENCE_realigned.vcf
# OUTPUT : nanopore_sv_calling/OAC_STRATFORD/OAC_STRATFORD_realigned.vcf
# OUTPUT : nanopore_sv_calling/OT09-03/OT09-03_realigned.vcf
# OUTPUT : nanopore_sv_calling/QS5091.50J/QS5091.50J_realigned.vcf
# OUTPUT : nanopore_sv_calling/ROLAND/ROLAND_realigned.vcf
			   output_vcf = paste0(i, "_realigned.vcf"),
			   ncores = 1,
			   reference_window = 500,
			   reads_window = 200,
			   min_overlap = 0.5,
			   min_identity = 85,
			   max_gaps = 15,
			   max_distance = 0.5,
			   max_offset = 50,
			   max_svlen = 50000,
			   # DEPENDENCY : scripts/breakpoint_refinement/age_realign.sh
			   age_script = "../../scripts/breakpoint_refinement/age_realign.sh",
			   samtools = samtools,
			   minimap2 = minimap2,
			   age = age,
			   wtdbg2 = wtdbg2,
			   wtpoa_cns = wtpoa_cns,
			   # DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
			   refgenome = "../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
# DEPENDENCY : nanopore_data/AC2001_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/ALTA_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/MAPLE_ISLE_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/MAPLE_PRESTO_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_09_35C_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_CARMAN_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_DRAYTON_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_EMBRO_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_LAKEVIEW_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_MADOC_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_OXFORD_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_PETREL_2_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_PRUDENCE_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OAC_STRATFORD_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/OT09-03_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/QS5091.50J_porechopped_aligned.sort.bam
# DEPENDENCY : nanopore_data/ROLAND_porechopped_aligned.sort.bam
			   nanopore_bam = paste0("../../nanopore_data/", i, "_porechopped_aligned.sort.bam"))

	setwd("..")
}

