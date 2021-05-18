
# Creating some variables for executables
# Using R 3.5.0 for figures because of a problem with resolution when using R 4.0
AGE = /home/malem420/programs/AGE/age_align
ASMVAR = /home/malem420/programs/AsmVar/src/AsmvarDetect/ASV_VariantDetector
BAYESTYPERTOOLS = /home/malem420/programs/bayesTyper_v1.5_linux_x86_64/bin/bayesTyperTools
BCFTOOLS = /home/malem420/programs/bcftools/bcftools
CIRCOS = /home/malem420/programs/circos-0.69-9/bin/circos
LASTAL = /home/malem420/programs/last-1047/src/lastal
LASTSPLIT = /home/malem420/programs/last-1047/src/last-split
MINIMAP2 = /home/malem420/programs/minimap2/minimap2
NGMLR = /home/malem420/programs/ngmlr/ngmlr-0.2.7/ngmlr
PORECHOP = /home/malem420/programs/Porechop/porechop-runner.py
R_FIG_COMMAND = /usr/bin/Rscript
R_RUN_COMMAND = /prg/R/4.0/bin/Rscript
SAMTOOLS = /home/malem420/programs/samtools/samtools
SOAPDENOVO2 = /prg/SOAPdenovo/2.04/SOAPdenovo-63mer
SNIFFLES = /home/malem420/programs/Sniffles-master/bin/sniffles-core-1.0.11/sniffles
WTDBG2 = /home/malem420/programs/wtdbg2/wtdbg2
WTPOA_CNS = /home/malem420/programs/wtdbg2/wtpoa-cns

# Creating some variables for more readable coding
# Supplemental figures S1 to S19 plus figures S20 and S21 which depend on other files
SUPFIGURES = $(shell seq 1 19 | xargs -I {} echo figures/figure_s{}.png) \
	     figures/Gm04_2257090_annotated.png \
	     figures/Gm04_2254504_annotated.png

SUPTABLES = $(shell seq 1 8 | xargs -I {} echo tables/table_s{}.csv)

FIGURES = $(shell seq 1 6 | xargs -I {} echo figures/figure_{}.png)

TABLES = $(shell seq 1 3 | xargs -I {} echo tables/table_{}.png)

# --- This target prepares all the figures, tables, and supplemental data
all: Supplemental_Data.pdf $(TABLES) $(FIGURES)

# --- The following section prepares de Supplemental Data file
Supplemental_Data.pdf : Supplemental_Data.tex references.bib genome_research.bst $(SUPFIGURES) $(SUPTABLES) 
	pdflatex Supplemental_Data.tex; bibtex Supplemental_Data; pdflatex Supplemental_Data.tex; pdflatex Supplemental_Data.tex

figures/figure_%.png : figures/figure_%.R
	cd figures; $(R_FIG_COMMAND) $(<F)

figures/figure_s1.png: figures/figure_s1.R \
	sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData \
	scripts/make_plot_data.R

tables/table_s1.csv tables/table_s2.csv tables/table_s3.csv: tables/formatting_sup_tables.R
	cd tables; $(R_RUN_COMMAND) formatting_sup_tables.R

tables/table_s4.csv: tables/table_s4.R
	cd tables; $(R_RUN_COMMAND) table_s4.R

tables/table_s5.csv: tables/table_s5.R
	cd tables; $(R_RUN_COMMAND) table_s5.R

tables/table_s6.csv tables/table_s7.csv tables/table_s8.csv: tables/tables_s6_s7_s8.R 
	cd tables; $(R_RUN_COMMAND) tables_s6_s7_s8.R

supfigures : $(SUPFIGURES)

suptables : $(SUPTABLES)

# --- The following section prepares the main tables for inclusion in the manuscript

tables: $(TABLES)

tables/table_1.png : tables/table_1.tex tables/table_1.csv
	cd tables; \
	pdflatex table_1.tex; \
	pdftoppm table_1.pdf -singlefile -r 500 -png table_1

tables/table_2.png : tables/table_2.tex tables/table_2.csv
	cd tables; \
	pdflatex table_2.tex; \
	pdftoppm table_2.pdf -singlefile -r 500 -png table_2

tables/table_3.png : tables/table_3.tex tables/table_3.csv
	cd tables; \
	pdflatex table_3.tex; \
	pdftoppm table_3.pdf -singlefile -r 500 -png table_3

tables/table_1.csv: tables/gather_table_1_data.sh
	cd tables; ./gather_table_1_data.sh

tables/table_2.csv: tables/table_2.R
	cd tables; $(R_RUN_COMMAND) table_2.R

tables/table_3.csv: tables/table_3.R
	cd tables; $(R_RUN_COMMAND) table_3.R

# --- The following section prepares the main figures for inclusion in the manuscript
figures : $(FIGURES)

figures/figure_1.png: figures/figure_1.R
	cd figures; $(R_FIG_COMMAND) figure_1.R

figures/figure_2.png: figures/figure_2.R
	cd figures; $(R_FIG_COMMAND) figure_2.R

# --- Section for the Circos figure
CIRCD = figures/figure_3_circos
EXTDIR = external/circos_config_files

figures/figure_3.png: $(CIRCD)/circos.conf $(EXTDIR)/housekeeping.conf $(EXTDIR)/housekeeping.conf $(CIRCD)/image.conf \
	$(CIRCD)/Gmax_karyotype.txt $(CIRCD)/dummy_karyotype.txt $(CIRCD)/ideogram.conf $(CIRCD)/ticks.conf $(CIRCD)/plots.conf \
	$(CIRCD)/highlights.conf $(CIRCD)/genes_heatmap.txt $(CIRCD)/snp_density.txt $(CIRCD)/sv_counts.txt $(CIRCD)/ref_ltr.txt \
	$(CIRCD)/poly_ltr.txt $(CIRCD)/ref_dna.txt $(CIRCD)/poly_dna.txt $(CIRCD)/legendA.txt $(CIRCD)/legendB.txt $(CIRCD)/legendC.txt \
	$(CIRCD)/legendD.txt $(CIRCD)/legendE.txt $(CIRCD)/dna_axis_min.txt $(CIRCD)/dna_axis_max.txt $(CIRCD)/ltr_axis_min.txt \
	$(CIRCD)/ltr_axis_max.txt $(CIRCD)/sv_highlights.txt $(CIRCD)/ltr_highlights.txt $(CIRCD)/dna_highlights.txt
	cd $(CIRCD); $(CIRCOS)

$(CIRCD)/Gmax_karyotype.txt: $(CIRCD)/make_karyotype_file.R refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd $(CIRCD); $(R_RUN_COMMAND) make_karyotype_file.R

$(CIRCD)/genes_heatmap.txt: $(CIRCD)/make_genes_heatmap_data.R $(CIRCD)/gmax4_3Mb_bins.RData
	cd $(CIRCD); $(R_RUN_COMMAND) make_genes_heatmap_data.R

$(CIRCD)/gmax4_3Mb_bins.RData: $(CIRCD)/make_3Mb_bins.R refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd $(CIRCD); $(R_RUN_COMMAND) make_3Mb_bins.R

$(CIRCD)/snp_density.txt: $(CIRCD)/make_snp_track.R structure_analysis/platypus_filtered_snps.vcf $(CIRCD)/gmax4_3Mb_bins.RData
	cd $(CIRCD); $(R_RUN_COMMAND) make_snp_track.R

$(CIRCD)/sv_counts.txt: $(CIRCD)/make_SV_track.R $(CIRCD)/gmax4_3Mb_bins.RData
	cd $(CIRCD); $(R_RUN_COMMAND) make_SV_track.R

$(CIRCD)/ref_ltr.txt $(CIRCD)/poly_ltr.txt $(CIRCD)/ref_dna.txt $(CIRCD)/poly_dna.txt: $(CIRCD)/make_te_tracks.R \
	refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3 te_analysis/polymorphic_tes.tsv $(CIRCD)/gmax4_3Mb_bins.RData
	cd $(CIRCD); $(R_RUN_COMMAND) make_te_tracks.R

$(CIRCD)/dna_axis_min.txt $(CIRCD)/dna_axis_max.txt $(CIRCD)/ltr_axis_min.txt $(CIRCD)/ltr_axis_max.txt: $(CIRCD)/make_axis_labels.R \
	$(CIRCD)/ref_ltr.txt $(CIRCD)/poly_ltr.txt $(CIRCD)/ref_dna.txt $(CIRCD)/poly_dna.txt
	cd $(CIRCD); $(R_RUN_COMMAND) make_axis_labels.R

$(CIRCD)/sv_highlights.txt $(CIRCD)/ltr_highlights.txt $(CIRCD)/dna_highlights.txt: $(CIRCD)/make_highlights_tracks.R \
	$(CIRCD)/sv_counts.txt $(CIRCD)/ref_ltr.txt $(CIRCD)/poly_ltr.txt $(CIRCD)/ref_dna.txt $(CIRCD)/poly_dna.txt
	cd $(CIRCD); $(R_RUN_COMMAND) make_highlights_tracks.R

# --- End of the Circos figure section

figures/figure_4.png: figures/figure_4.R
	cd figures; $(R_FIG_COMMAND) figure_4.R

figures/figure_5.png: figures/figure_5.R
	cd figures; $(R_FIG_COMMAND) figure_5.R

figures/figure_6.png: figures/figure_6.R
	cd figures; $(R_FIG_COMMAND) figure_6.R

# --- This section processes the raw basecalled Nanopore reads using Porechop and aligns them to the reference genome
NANOPORE_READS := $(shell cat utilities/flowcell_names.txt | xargs -I {} echo nanopore_data/{}.fastq.gz)

nanopore_data/NANOPORE_ALIGNMENT : $(NANOPORE_READS) \
	nanopore_data/read_processing_alignment.sh \
	utilities/flowcell_names.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd nanopore_data ; ./read_processing_alignment.sh $(PORECHOP) $(NGMLR) $(SAMTOOLS) ; touch NANOPORE_ALIGNMENT

# --- This section prepares the Oxford Nanopore Sniffles VCF files from the sorted bam files
NANOPORE_SORTED_BAM := $(shell tail -n+2 utilities/line_ids.txt | cut -f1 | xargs -I {} echo nanopore_data/{}_porechopped_aligned.sort.bam)
NANOPORE_FILTERED_SVS := $(shell tail -n+2 utilities/line_ids.txt | cut -f1 | xargs -I {} echo nanopore_sv_calling/{}_hom70_filtered.vcf)
NANOPORE_REFINED_SVS := $(shell tail -n+2 utilities/line_ids.txt | cut -f1 | xargs -I {} echo nanopore_sv_calling/{}/{}_realigned.vcf)
NANOPORE_NORMALIZED_SVS := $(shell tail -n+2 utilities/line_ids.txt | cut -f1 | xargs -I {} echo nanopore_sv_calling/{}_normalized_ids.vcf)

# Calling SVs with Sniffles and filtering the output
nanopore_sv_calling/SV_CALLING : nanopore_data/NANOPORE_ALIGNMENT $(NANOPORE_SORTED_BAM) \
	nanopore_sv_calling/sniffles_calling.sh \
	nanopore_sv_calling/filter_sniffles_vcfs.R \
	scripts/filter_sniffles.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd nanopore_sv_calling ; ./sniffles_calling.sh $(SNIFFLES) $(R_RUN_COMMAND) ; touch SV_CALLING

# Refining the SV breakpoints
nanopore_sv_calling/SV_REFINEMENT : nanopore_sv_calling/SV_CALLING $(NANOPORE_FILTERED_SVS) \
	nanopore_data/NANOPORE_ALIGNMENT $(NANOPORE_SORTED_BAM) \
	scripts/breakpoint_refinement/refine_breakpoints.R \
	scripts/breakpoint_refinement/gather_align_data.R \
	scripts/breakpoint_refinement/parse_age.R \
	scripts/breakpoint_refinement/parse_svinfo.R \
	scripts/breakpoint_refinement/revcomp.R \
	scripts/breakpoint_refinement/call_age.R \
	scripts/breakpoint_refinement/update_breakpoints.R \
	utilities/line_ids.txt \
	scripts/breakpoint_refinement/age_realign.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd nanopore_sv_calling ; $(R_RUN_COMMAND) refine_breakpoints.R $(SAMTOOLS) $(MINIMAP2) $(AGE) $(WTDBG2) $(WTPOA_CNS) ; touch SV_REFINEMENT

# Processing the refined VCFs to prepare them for input to SVmerge;
# These files are also the reference Oxford Nanpore SVs for benchmarking
nanopore_sv_calling/SV_NORMALIZATION : nanopore_sv_calling/SV_REFINEMENT $(NANOPORE_REFINED_SVS) \
	nanopore_sv_calling/process_vcf_files.sh \
	nanopore_sv_calling/add_metainfo_all.R \
	nanopore_sv_calling/fix_vcfs.R \
	scripts/add_metainfo.R \
	scripts/fix_sniffles.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd nanopore_sv_calling ; ./process_vcf_files.sh $(BCFTOOLS) $(R_RUN_COMMAND) ; touch SV_NORMALIZATION

# --- This section calls and filters the SVs from AsmVar following de novo assembly with SOAPdenovo2
ILLUMINA_SINGLE_TRIMMED = $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_data/trimmed_fastq/{}/{}_sing_trimmed.fastq.gz)
FLASH_MERGED = $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/flash_merging/{}/out.extendedFrags.fastq.gz)
FLASH_UNMERGED = $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/flash_merging/{}/out.notCombined_1.fastq.gz)
LAST_ALIGNMENTS := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/last_alignment/{}_soapdenovo2_fm.maf)
SOAPDENOVO_ASSEMBLIES :=  $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/soap_assembly/{}/{}.contig)
ASMVAR_VCFS := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/asmvar_calling/{}/asmvar_results_Gm01.vcf)

# Assembling the sequences using SOAPdenovo2
illumina_sv_calling/asmvar/SOAPDENOVO_ASSEMBLY : $(ILLUMINA_SINGLE_TRIMMED) $(FLASH_SEQUENCES) \
	illumina_sv_calling/asmvar/assembly.sh \
	illumina_sv_calling/asmvar/soap_assembly/make_config.sh \
	illumina_sv_calling/asmvar/soap_assembly/config_template.txt \
	utilities/all_lines.txt
	cd illumina_sv_calling/asmvar ; ./assembly.sh $(SOAPDENOVO2) ; touch SOAPDENOVO_ASSEMBLY

# Aligning the de novo assemblies to the reference genome using LAST
illumina_sv_calling/asmvar/LAST_ALIGNMENT : illumina_sv_calling/asmvar/SOAPDENOVO_ASSEMBLY $(SOAPDENOVO_ASSEMBLIES) \
	illumina_sv_calling/asmvar/last_alignment.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd illumina_sv_calling/asmvar ; ./last_alignment.sh $(LASTAL) $(LASTSPLIT) ; touch LAST_ALIGNMENT

# Calling SVs with AsmVar on the aligned assemblies
illumina_sv_calling/asmvar/ASMVAR_CALLING : illumina_sv_calling/asmvar/LAST_ALIGNMENT $(LAST_ALIGNMENTS) \
	$(SOAPDENOVO_ASSEMBLIES) \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd illumina_sv_calling/asmvar ; ./asmvar_call.sh $(ASMVAR) ; touch ASMVAR_CALLING

# Filtering the VCFs resulting from calling AsmVar on assemblies
illumina_sv_calling/asmvar/asmvar_filtering/asmvar_svs.vcf : illumina_sv_calling/asmvar/ASMVAR_CALLING $(ASMVAR_VCFS) \
	illumina_sv_calling/asmvar/asmvar_filter.sh \
	utilities/all_lines.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai \
	scripts/add_svtype.awk \
	illumina_sv_calling/asmvar/svtype_header_line.txt \
	scripts/extract_svs_50.awk
	cd illumina_sv_calling/asmvar ; ./asmvar_filter.sh $(BCFTOOLS) $(BAYESTYPERTOOLS)

# --- The next section prepares the Illumina SV benchmarks from the Paragraph vcfs
ILLUMINA_BENCHMARK_VCFS := $(shell tail -n+2 utilities/line_ids.txt | cut -f2 | xargs -I {} echo sv_genotyping/illumina_svs/{}_results/genotypes.vcf.gz)
# Benchmark of Illumina SVs in non-repeat regions
sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData: \
	nanopore_sv_calling/SV_NORMALIZATION \
	sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_benchmark.R \
	utilities/line_ids.txt \
	$(ILLUMINA_BENCHMARK_VCFS) \
	$(NANOPORE_NORMALIZED_SVS) \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/repeat_regions/non_repeated_regions.bed \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; $(R_RUN_COMMAND) norepeat_benchmark.R


