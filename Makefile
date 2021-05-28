# Creating some variables for executables
# Using R 3.5.0 for figures because of a problem with resolution when using R 4.0
AGE = /home/malem420/programs/AGE/age_align
ASMVAR = /home/malem420/programs/AsmVar/src/AsmvarDetect/ASV_VariantDetector
BAYESTYPER = /home/malem420/programs/bayesTyper_v1.5_linux_x86_64/bin/bayesTyper
BAYESTYPERTOOLS = /home/malem420/programs/bayesTyper_v1.5_linux_x86_64/bin/bayesTyperTools
BBDUK = /home/malem420/programs/bbmap/bbduk.sh
BAMADDRG = /prg/bamaddrg/1.0/bamaddrg
BGZIP = /prg/htslib/1.10.2/bin/bgzip
BWA = /prg/bwa/0.7.17/bwa
BCFTOOLS = /home/malem420/programs/bcftools/bcftools
BCFTOOLS_PLUGIN_PATH = /home/malem420/programs/bcftools/plugins
BLASTN = /home/malem420/programs/ncbi-blast-2.11.0+/bin/blastn
CIRCOS = /home/malem420/programs/circos-0.69-9/bin/circos
FASTSTRUCTURE = /prg/fastStructure/1.0/structure.py
FLASH = /home/malem420/programs/FLASH-1.2.11-Linux-x86_64/flash
IDXDEPTH = /home/malem420/programs/paragraph/bin/idxdepth
KMC = /home/malem420/programs/KMC3.linux/kmc
LASTAL = /home/malem420/programs/last-1047/src/lastal
LASTSPLIT = /home/malem420/programs/last-1047/src/last-split
MAKEBLASTDB = /home/malem420/programs/ncbi-blast-2.11.0+/bin/makeblastdb
MANTA = /home/malem420/programs/manta-1.6.0.centos6_x86_64/bin/configManta.py
MINIMAP2 = /home/malem420/programs/minimap2/minimap2
MULTIGRMPY = /home/malem420/programs/paragraph/bin/multigrmpy.py
NANOPLOT = /home/malem420/.local/bin/NanoPlot
NGMLR = /home/malem420/programs/ngmlr/ngmlr-0.2.7/ngmlr
PLATYPUS = /prg/platypus/0.8.1.1/bin/Platypus.py
PLINK = /prg/plink/1.90b5.3/bin/plink
PORECHOP = /home/malem420/programs/Porechop/porechop-runner.py
R_FIG_COMMAND = /usr/bin/Rscript
R_RUN_COMMAND = /prg/R/4.0/bin/Rscript
SAMTOOLS = /home/malem420/programs/samtools/samtools
SMOOVE = /home/malem420/programs/smoove/smoove
SOAPDENOVO2 = /prg/SOAPdenovo/2.04/SOAPdenovo-63mer
SNIFFLES = /home/malem420/programs/Sniffles-master/bin/sniffles-core-1.0.11/sniffles
SVABA = /home/malem420/programs/svaba/bin/svaba
SVMERGE = /home/malem420/programs/SVanalyzer-install/bin/SVmerge
TABIX = /prg/htslib/1.10.2/bin/tabix
VCFTOOLS = /prg/vcftools/0.1.16/bin/vcftools
VG = /home/malem420/programs/vg/vg
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

figures/figure_s1.png : sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData scripts/make_plot_data.R

figures/figure_s2.png : sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData scripts/make_plot_data.R

figures/figure_s3.png : sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData scripts/make_plot_data.R

figures/figure_s4.png : depth_distributions/average_depth.RData \
	utilities/line_ids.txt \
	sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData

figures/figure_s5.png : sv_genotyping/illumina_svs/sveval_benchmarks/frequency_RData/sveval_frequency_rates.RData scripts/make_plot_data.R

figures/figure_s6.png : sv_genotyping/illumina_svs/sveval_benchmarks/frequency_RData/sveval_frequency_rates.RData scripts/make_plot_data.R

figures/figure_s7.png : sv_genotyping/illumina_svs/sveval_benchmarks/NCALLERS_ILLUMINA_BENCHMARK scripts/make_plot_data.R

figures/figure_s8.png : sv_genotyping/illumina_svs/size_distribution.tsv

figures/figure_s9.png : sv_genotyping/illumina_svs/sveval_benchmarks/NCALLERS_ILLUMINA_BENCHMARK scripts/format_sveval_plotting_data.R

figures/figure_s10.png : sv_genotyping/illumina_svs/sveval_benchmarks/NCALLERS_ILLUMINA_BENCHMARK scripts/format_sveval_plotting_data.R

figures/figure_s11.png : nanoplot_stats/N50_stats.txt \
	utilities/line_ids.txt \
	sv_genotyping/nanopore_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData

figures/figure_s12.png : sv_genotyping/nanopore_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData scripts/make_plot_data.R

figures/figure_s13.png : sv_genotyping/nanopore_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData scripts/make_plot_data.R

figures/figure_s14.png : sv_genotyping/nanopore_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData scripts/make_plot_data.R

figures/figure_s15.png : sv_genotyping/combined_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData scripts/make_plot_data.R

figures/figure_s16.png : sv_genotyping/combined_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData scripts/make_plot_data.R

figures/figure_s17.png : breakpoint_refinement_analysis/raw_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData \
	breakpoint_refinement_analysis/refined_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData \
	scripts/make_plot_data.R

figures/figure_s18.png : breakpoint_refinement_analysis/raw_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData \
	breakpoint_refinement_analysis/refined_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData \
	scripts/make_plot_data.R

figures/figure_s19.png : breakpoint_refinement_analysis/raw_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData \
	breakpoint_refinement_analysis/refined_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData \
	scripts/make_plot_data.R

tables/table_%.csv : tables/table_%.R
	cd tables ; $(R_RUN_COMMAND) $(<F)

tables/table_s1.csv : tables/lab_methods_table.csv

tables/table_s2.csv : tables/lab_methods_table.csv

tables/table_s3.csv : tables/table_s3_data.txt

tables/table_s4.csv : nanopore_sv_calling/all_metainfo.RData

tables/table_s5.csv : gene_analysis/allele_frequency_permutations.RData

tables/table_s6.csv : gene_analysis/GO_ANALYSIS scripts/format_go_table.R

tables/table_s7.csv : gene_analysis/GO_ANALYSIS scripts/format_go_table.R

tables/table_s8.csv : gene_analysis/GO_ANALYSIS scripts/format_go_table.R

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

figures/figure_1.png : sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData scripts/make_plot_data.R

figures/figure_2.png : sv_genotyping/nanopore_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData scripts/make_plot_data.R 

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

$(CIRCD)/genes_heatmap.txt: $(CIRCD)/make_genes_heatmap_data.R $(CIRCD)/gmax4_3Mb_bins.RData gene_analysis/GENE_OVERLAP_ANALYSIS
	cd $(CIRCD); $(R_RUN_COMMAND) make_genes_heatmap_data.R

$(CIRCD)/gmax4_3Mb_bins.RData: $(CIRCD)/make_3Mb_bins.R refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd $(CIRCD); $(R_RUN_COMMAND) make_3Mb_bins.R

$(CIRCD)/snp_density.txt: $(CIRCD)/make_snp_track.R structure_analysis/platypus_filtered_snps.vcf $(CIRCD)/gmax4_3Mb_bins.RData
	cd $(CIRCD); $(R_RUN_COMMAND) make_snp_track.R

$(CIRCD)/sv_counts.txt: $(CIRCD)/make_SV_track.R $(CIRCD)/gmax4_3Mb_bins.RData sv_genotyping/combined_svs/combined_paragraph_filtered.vcf
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

# --- This section creates the bed file of repeats from the reference genome and Phytozome repeat annotation
refgenome/repeat_regions/non_repeated_regions.bed : refgenome/repeat_regions/make_repeat_bed.R \
	refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3 \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd refgenome/repeat_regions ; $(R_RUN_COMMAND) make_repeat_bed.R

# --- This section compute the average sequencing depths used for plotting figure S4
depth_distributions/average_depth.RData : \
	nanopore_data/NANOPORE_ALIGNMENT $(NANOPORE_SORTED_BAM) \
	utilities/line_ids.txt \
	depth_distributions/depth_distributions.sh \
	scripts/depth_distribution.sh \
	depth_distributions/compute_average_depth.R
	cd depth_distributions/ ; ./depth_distributions.sh $(SAMTOOLS) ; $(R_RUN_COMMAND) compute_average_depth.R

# --- Generating the SV size distributions necessary for figure S8
sv_genotyping/illumina_svs/size_distribution.tsv : sv_genotyping/illumina_svs/svmerged.clustered.vcf \
	sv_genotyping/illumina_svs/extract_size_distribution.sh
	cd sv_genotyping/illumina_svs/ ; ./extract_size_distribution.sh $(BCFTOOLS)

# --- Generating the Nanoplot N50 data for figure S11
nanoplot_stats/NANOPLOT_BAM_STATS : nanopore_data/NANOPORE_ALIGNMENT $(NANOPORE_SORTED_BAM) \
	nanoplot_stats/nanoplot_sorted_bam_all.sh \
	utilities/line_ids.txt
	cd nanoplot_stats ; ./nanoplot_sorted_bam_all.sh $(NANOPLOT) ; touch nanoplot_stats/NANOPLOT_BAM_STATS

nanoplot_stats/N50_stats.txt : nanoplot_stats/NANOPLOT_BAM_STATS \
	nanoplot_stats/gather_N50_stats.sh \
	utilities/line_ids.txt
	cd nanoplot_stats ; ./gather_N50_stats.sh > N50_stats.txt

# Generating the data for Table S3
tables/table_s3_data.txt : utilities/line_ids.txt \
       	tables/gather_table_s3_data.sh \
       	nanopore_sv_calling/SV_NORMALIZATION \
	scripts/count_svtypes_svsizes.awk
	cd tables ; ./gather_table_s3_data.sh

# Generation the data for Table S4
nanopore_sv_calling/all_metainfo.RData : nanopore_sv_calling/SV_NORMALIZATION \
	nanopore_sv_calling/gather_metainfo.R
	cd nanopore_sv_calling/ ; $(R_RUN_COMMAND) gather_metainfo.R $(BCFTOOLS)

# Generating the data for Table S5
gene_analysis/allele_frequency_permutations.RData : gene_analysis/allele_frequency_permutations.R \
	gene_analysis/GENE_OVERLAP_ANALYSIS
	cd gene_analysis ; $(R_RUN_COMMAND) allele_frequency_permutations.R

# GENE_OVERLAP_ANALYSIS includes :
# gene_analysis/overlap_data.txt
# gene_analysis/overlap_norepeat.txt
# gene_analysis/genes.RData
# gene_analysis/cds.RData
# gene_analysis/upstream5kb.RData
# gene_analysis/feature_widths.RData
# gene_analysis/norepeat_widths.RData
gene_analysis/GENE_OVERLAP_ANALYSIS : sv_genotyping/combined_svs/combined_paragraph_filtered.vcf \
	gene_analysis/gene_overlap_analysis.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/Gmax_508_Wm82.a4.v1.gene_exons.gff3 \
	refgenome/repeat_regions/non_repeated_regions.bed
	cd gene_analysis/ ; $(R_RUN_COMMAND) gene_overlap_analysis.R

# Generating the data for tables S6, S7 and S8
# Outputs the following files :
# bp_over_summary.RData
# bp_under_summary.RData
# pfam_over_summary.RData
gene_analysis/GO_ANALYSIS : gene_analysis/GENE_OVERLAP_ANALYSIS \
	gene_analysis/go_analysis.R \
	gene_analysis/soybase_genome_annotation_v4.0_04-20-2021.txt
	cd gene_analysis/ ; $(R_RUN_COMMAND) go_analysis.R

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

# --- This section takes care of Illumina read trimming and alignment
ILLUMINA_READS_1 = $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_data/raw_fastq/{}_R1.fastq.gz)
ILLUMINA_READS_2 = $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_data/raw_fastq/{}_R2.fastq.gz)
ILLUMINA_PAIRED_TRIMMED := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_data/trimmed_fastq/{}/{}_R1_trimmed.fastq.gz)
ILLUMINA_SINGLE_TRIMMED := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_data/trimmed_fastq/{}/{}_sing_trimmed.fastq.gz)
ILLUMINA_ALIGNED_READS := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_data/aligned_reads/{}/{}_all.sort.bam)

illumina_data/ILLUMINA_TRIMMING : $(ILLUMINA_READS_1) $(ILLUMINA_READS_2) \
	illumina_data/bbduk_trimming.sh \
	utilities/all_lines.txt
	cd illumina_data ; ./bbduk_trimming.sh $(BBDUK) ; touch illumina_data/ILLUMINA_TRIMMING

illumina_data/ILLUMINA_ALIGNMENT : illumina_data/ILLUMINA_TRIMMING $(ILLUMINA_PAIRED_TRIMMED) $(ILLUMINA_SINGLE_TRIMMED) \
	illumina_data/bwa_alignment.sh \
	utilities/all_lines.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd illumina_data ; ./bwa_alignment.sh $(BWA) $(SAMTOOLS) $(BAMADDRG) ; touch ILLUMINA_ALIGNMENT

# --- This section calls and filters the SVs from AsmVar following de novo assembly with SOAPdenovo2
FLASH_MERGED := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/flash_merging/{}/out.extendedFrags.fastq.gz)
FLASH_UNMERGED := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/flash_merging/{}/out.notCombined_1.fastq.gz)
LAST_ALIGNMENTS := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/last_alignment/{}_soapdenovo2_fm.maf)
SOAPDENOVO_ASSEMBLIES :=  $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/soap_assembly/{}/{}.contig)
ASMVAR_VCFS := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/asmvar/asmvar_calling/{}/asmvar_results_Gm01.vcf)

# Merging the trimmed sequences with FLASH
illumina_sv_calling/asmvar/FLASH_MERGING : illumina_data/ILLUMINA_TRIMMING $(ILLUMINA_PAIRED_TRIMMED) \
	illumina_sv_calling/asmvar/flash_merging.sh
	cd illumina_sv_calling/asmvar ; ./flash_merging.sh $(FLASH) ; touch FLASH_MERGING

# Assembling the sequences using SOAPdenovo2
illumina_sv_calling/asmvar/SOAPDENOVO_ASSEMBLY : illumina_data/ILLUMINA_TRIMMING $(ILLUMINA_SINGLE_TRIMMED) \
	illumina_sv_calling/asmvar/FLASH_MERGING $(FLASH_SEQUENCES) \
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
	illumina_sv_calling/asmvar/SOAPDENOVO_ASSEMBLY $(SOAPDENOVO_ASSEMBLIES) \
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

# --- This section calls and filters the SVs discovered using manta
MANTA_CANDIDATE_SVS := $(shell seq 1 10 | xargs -I {} echo illumina_sv_calling/manta/manta_sample{}/MantaWorkflow/results/variants/candidateSV.vcf.gz)

illumina_sv_calling/manta/MANTA_CALLING : illumina_data/ILLUMINA_ALIGNMENT $(ILLUMINA_ALIGNED_READS) \
	illumina_sv_calling/manta/manta_call.sh \
	illumina_sv_calling/manta/manta_config.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd illumina_sv_calling/manta ; ./manta_call.sh $(MANTA) ; touch MANTA_CALLING

illumina_sv_calling/manta/manta_svs.vcf : illumina_sv_calling/manta/MANTA_CALLING $(MANTA_CANDIDATE_SVS) \
	illumina_sv_calling/manta/manta_filter.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	illumina_sv_calling/manta/ACO_header_line.txt \
	scripts/extract_svs_50.awk
	cd illumina_sv_calling/manta ; ./manta_filter.sh $(BCFTOOLS) $(BAYESTYPERTOOLS) $(TABIX)

# --- This section calls and filters the SVs discovered using SvABA
SVABA_INDEL_VCFS := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/svaba/{}.svaba.indel.vcf)
SVABA_SV_VCFS := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_sv_calling/svaba/{}.svaba.sv.vcf)

illumina_sv_calling/svaba/SVABA_CALLING : illumina_data/ILLUMINA_ALIGNMENT $(ILLUMINA_ALIGNED_READS) \
	illumina_sv_calling/svaba/svaba_call.sh \
	utilities/all_lines.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd illumina_sv_calling/svaba ; ./svaba_call.sh $(SVABA) ; touch SVABA_CALLING

illumina_sv_calling/svaba/svaba_svs.vcf : illumina_sv_calling/svaba/SVABA_CALLING $(SVABA_INDEL_VCFS) $(SVABA_SV_VCFS) \
	illumina_sv_calling/svaba/convert_svaba.R \
	illumina_sv_calling/svaba/svaba_filter.sh \
	scripts/svaba_process.R \
	utilities/all_lines.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai \
	illumina_sv_calling/svaba/annotate_svtype.awk \
	illumina_sv_calling/svaba/ACO_header_line.txt \
	scripts/extract_svs_50.awk
	cd illumina_sv_calling/svaba ; $(R_RUN_COMMAND) convert_svaba.R ; ./svaba_filter.sh $(BGZIP) $(TABIX) $(BCFTOOLS)

# --- This section calls and filters the SVs discovered using smoove
illumina_sv_calling/smoove/all_samples.smoove.square.vcf.gz : illumina_data/ILLUMINA_ALIGNMENT $(ILLUMINA_ALIGNED_READS) \
	illumina_sv_calling/smoove/smoove_call.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd illumina_sv_calling/smoove ; ./smoove_call.sh $(SMOOVE) $(TABIX)

illumina_sv_calling/smoove/smoove_svs.vcf : illumina_sv_calling/smoove/all_samples.smoove.square.vcf.gz \
	illumina_sv_calling/smoove/smoove_filter.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	illumina_sv_calling/smoove/ACO_header_line.txt \
	scripts/extract_svs_50.awk
	cd illumina_sv_calling/smoove ; ./smoove_filter.sh $(BCFTOOLS) $(BAYESTYPERTOOLS)

# --- This section prepares the files that are needed for all Paragraph runs, irrespective of the SV dataset
PARAGRAPH_MANIFEST_FILES := $(shell cat utilities/all_lines.txt | xargs -I {} echo sv_genotyping/manifest_files/{}_manifest.txt)

sv_genotyping/MANIFEST_FILES : illumina_data/ILLUMINA_ALIGNMENT $(ILLUMINA_ALIGNED_READS) \
	sv_genotyping/prepare_manifest_files.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt \
	utilities/read_lengths.txt
	cd sv_genotyping ; ./prepare_manifest_files.sh $(IDXDEPTH) ; touch MANIFEST_FILES

# --- This section genotypes the Illumina SVs using Paragraph, first preparing them with SVmerge
sv_genotyping/illumina_svs/svmerged.clustered.vcf : sv_genotyping/illumina_svs/svmerge_files.txt \
	illumina_sv_calling/asmvar/asmvar_filtering/asmvar_svs.vcf \
	illumina_sv_calling/manta/manta_svs.vcf \
	illumina_sv_calling/smoove/smoove_svs.vcf \
	illumina_sv_calling/svaba/svaba_svs.vcf \
	sv_genotyping/illumina_svs/SVmerge.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd sv_genotyping/illumina_svs ; ./SVmerge.sh $(SVMERGE)

sv_genotyping/illumina_svs/PARAGRAPH_ILLUMINA_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT $(ILLUMINA_ALIGNED_READS) \
	sv_genotyping/MANIFEST_FILES $(PARAGRAPH_MANIFEST_FILES) \
	sv_genotyping/illumina_svs/svmerged.clustered.vcf \
	sv_genotyping/illumina_svs/run_paragraph.sh \
	scripts/addMissingPaddingGmax4.py \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd sv_genotyping/illumina_svs ; ./run_paragraph.sh $(BCFTOOLS) $(BGZIP) $(TABIX) $(MULTIGRMPY) ; touch PARAGRAPH_ILLUMINA_GENOTYPING

# --- This section genotypes the Oxford Nanopore SVs using Paragraph, first preparing them with SVmerge
sv_genotyping/nanopore_svs/svmerged_clustered_sorted.vcf : sv_genotyping/nanopore_svs/svmerge_files.txt \
	nanopore_sv_calling/SV_NORMALIZATION $(NANOPORE_NORMALIZED_SVS) \
	sv_genotyping/nanopore_svs/SVmerge.sh \
	sv_genotyping/nanopore_svs/merge_realigned.R \
	scripts/merge_realigned_variants.R \
	sv_genotyping/nanopore_svs/svmerged_sort.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd sv_genotyping/nanopore_svs ; ./SVmerge.sh $(SVMERGE) ; $(R_RUN_COMMAND) merge_realigned.R ; ./svmerged_sort.sh $(BCFTOOLS)

sv_genotyping/nanopore_svs/PARAGRAPH_NANOPORE_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT $(ILLUMINA_ALIGNED_READS) \
	sv_genotyping/MANIFEST_FILES $(PARAGRAPH_MANIFEST_FILES) \
	sv_genotyping/nanopore_svs/svmerged_clustered_sorted.vcf \
	sv_genotyping/nanopore_svs/run_paragraph.sh \
	scripts/addMissingPaddingGmax4.py \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd sv_genotyping/nanopore_svs ; ./run_paragraph.sh $(BCFTOOLS) $(BGZIP) $(TABIX) $(MULTIGRMPY) ; touch PARAGRAPH_NANOPORE_GENOTYPING

# --- This section genotypes the combined Illumina/Oxford Nanopore SVs using Paragraph, first preparing them with SVmerge
sv_genotyping/combined_svs/illumina_merged_sorted.vcf : sv_genotyping/combined_svs/svmerge_files.txt \
	sv_genotyping/illumina_svs/svmerged.clustered.vcf \
	sv_genotyping/nanopore_svs/svmerged_clustered_sorted.vcf \
	sv_genotyping/combined_svs/SVmerge.sh \
	sv_genotyping/combined_svs/header_lines.txt \
	sv_genotyping/combined_svs/select_svs.R \
	sv_genotyping/combined_svs/sort_vcfs.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta 
	cd sv_genotyping/combined_svs ; ./SVmerge.sh $(BCFTOOLS) $(SVMERGE) ; $(R_RUN_COMMAND) select_svs.R ; ./sort_vcfs.sh $(BCFTOOLS)

sv_genotyping/combined_svs/PARAGRAPH_COMBINED_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT $(ILLUMINA_ALIGNED_READS) \
	sv_genotyping/MANIFEST_FILES $(PARAGRAPH_MANIFEST_FILES) \
	sv_genotyping/combined_svs/illumina_merged_sorted.vcf \
	sv_genotyping/combined_svs/run_paragraph.sh \
	scripts/addMissingPaddingGmax4.py \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd sv_genotyping/combined_svs ; ./run_paragraph.sh $(BCFTOOLS) $(BGZIP) $(TABIX) $(MULTIGRMPY) ; touch PARAGRAPH_COMBINED_GENOTYPING

# Merging the VCF files of the combined Illumina/Nanopore SVs genotyped by Paragraph
sv_genotyping/combined_svs/combined_paragraph_merged.vcf : sv_genotyping/combined_svs/PARAGRAPH_COMBINED_GENOTYPING \
	sv_genotyping/combined_svs/merge_vcf_files.sh \
	utilities/all_lines.txt
	cd sv_genotyping/combined_svs ; ./merge_vcf_files.sh $(BCFTOOLS) $(TABIX)

# Filtering the merged vcf file for minimum number of supporting reads, minimum homozygous ALT count, and removing inversions/duplications
sv_genotyping/combined_svs/combined_paragraph_filtered.vcf : sv_genotyping/combined_svs/combined_paragraph_merged.vcf \
	sv_genotyping/combined_svs/filter_merged_vcf.sh
	cd sv_genotyping/combined_svs ; ./filter_merged_vcf.sh $(BCFTOOLS) $(BCFTOOLS_PLUGIN_PATH)

# --- This section prepares the Illumina SV benchmarks from the Paragraph vcfs
PARAGRAPH_ILLUMINA_VCFS := $(shell cat utilities/all_lines.txt | xargs -I {} echo sv_genotyping/illumina_svs/{}_results/genotypes.vcf.gz)

# Benchmark of Illumina SVs genome-wide
sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData : \
	sv_genotyping/illumina_svs/PARAGRAPH_ILLUMINA_GENOTYPING $(PARAGRAPH_ILLUMINA_VCFS) \
	nanopore_sv_calling/SV_NORMALIZATION $(NANOPORE_NORMALIZED_SVS) \
	sv_genotyping/illumina_svs/sveval_benchmarks/benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; $(R_RUN_COMMAND) benchmark.R

# Benchmark of Illumina SVs in non-repeat regions
sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData : \
	sv_genotyping/illumina_svs/PARAGRAPH_ILLUMINA_GENOTYPING $(PARAGRAPH_ILLUMINA_VCFS) \
	nanopore_sv_calling/SV_NORMALIZATION $(NANOPORE_NORMALIZED_SVS) \
	sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/repeat_regions/non_repeated_regions.bed \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; $(R_RUN_COMMAND) norepeat_benchmark.R

# Benchmark of the Illumina SVs as a function of the homozygous ALT count
# For this we first need to prepare the input vcf file for the benchmark
sv_genotyping/illumina_svs/illumina_paragraph_merged.vcf : sv_genotyping/illumina_svs/PARAGRAPH_ILLUMINA_GENOTYPING $(PARAGRAPH_ILLUMINA_VCFS) \
	sv_genotyping/illumina_svs/merge_vcf_files.sh \
	utilities/all_lines.txt
	cd sv_genotyping/illumina_svs ; ./merge_vcf_files.sh $(BCFTOOLS) $(TABIX)

sv_genotyping/illumina_svs/sveval_benchmarks/paragraph_svs_minDP2_annotated.vcf : \
	sv_genotyping/illumina_svs/illumina_paragraph_merged.vcf \
	sv_genotyping/illumina_svs/sveval_benchmarks/annotate_filter_population_svs.sh
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; ./annotate_filter_population_svs.sh $(BCFTOOLS) $(BCFTOOLS_PLUGIN_PATH)

# Computing the benchmarks using the number of ALT alleles in homozygous genotype calls (homozygous ALT counts) for the precision-recall curves
sv_genotyping/illumina_svs/sveval_benchmarks/frequency_RData/sveval_frequency_rates.RData : \
	sv_genotyping/illumina_svs/sveval_benchmarks/paragraph_svs_minDP2_annotated.vcf \
	nanopore_sv_calling/SV_NORMALIZATION $(NANOPORE_NORMALIZED_SVS) \
	sv_genotyping/illumina_svs/sveval_benchmarks/allele_frequency_benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; $(R_RUN_COMMAND) allele_frequency_benchmark.R

# Computing the benchmarks using the number of SV calling tools for the precision-recall curves
sv_genotyping/illumina_svs/sveval_benchmarks/paragraph_svs_minAC4_ncallers.vcf : \
	sv_genotyping/illumina_svs/sveval_benchmarks/paragraph_svs_minDP2_annotated.vcf \
	sv_genotyping/illumina_svs/sveval_benchmarks/filter_minAC_annotateCallers.sh \
	sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_header_line.txt
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; ./filter_minAC_annotateCallers.sh $(BCFTOOLS)

sv_genotyping/illumina_svs/sveval_benchmarks/NCALLERS_ILLUMINA_BENCHMARK : \
	sv_genotyping/illumina_svs/sveval_benchmarks/paragraph_svs_minAC4_ncallers.vcf \
	nanopore_sv_calling/SV_NORMALIZATION $(NANOPORE_NORMALIZED_SVS) \
	sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; $(R_RUN_COMMAND) ncallers_benchmark.R

# --- The next section prepares the Oxford Nanopore SV benchmarks from the Paragraph vcfs
PARAGRAPH_NANOPORE_VCFS := $(shell cat utilities/all_lines.txt | xargs -I {} echo sv_genotyping/nanopore_svs/{}_results/genotypes.vcf.gz)

# Benchmark of Oxford Nanopore SVs genome-wide
sv_genotyping/nanopore_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData : \
	sv_genotyping/nanopore_svs/PARAGRAPH_NANOPORE_GENOTYPING $(PARAGRAPH_NANOPORE_VCFS) \
	nanopore_sv_calling/SV_NORMALIZATION $(NANOPORE_NORMALIZED_SVS) \
	sv_genotyping/nanopore_svs/sveval_benchmarks/benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/nanopore_svs/sveval_benchmarks ; $(R_RUN_COMMAND) benchmark.R

# Benchmark of Oxford Nanopore SVs in non-repeat regions
sv_genotyping/nanopore_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData : \
	sv_genotyping/nanopore_svs/PARAGRAPH_NANOPORE_GENOTYPING $(PARAGRAPH_NANOPORE_VCFS) \
	nanopore_sv_calling/SV_NORMALIZATION $(NANOPORE_NORMALIZED_SVS) \
	sv_genotyping/nanopore_svs/sveval_benchmarks/norepeat_benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/repeat_regions/non_repeated_regions.bed \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/nanopore_svs/sveval_benchmarks ; $(R_RUN_COMMAND) norepeat_benchmark.R

# --- The next section prepares the combined Illumina/Oxford Nanopore SV benchmarks from the Paragraph vcfs
PARAGRAPH_COMBINED_VCFS := $(shell cat utilities/all_lines.txt | xargs -I {} echo sv_genotyping/combined_svs/{}_results/genotypes.vcf.gz)

# Benchmark of combined Illumina/Oxford Nanopore SVs genome-wide
sv_genotyping/combined_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData : \
	sv_genotyping/combined_svs/PARAGRAPH_COMBINED_GENOTYPING $(PARAGRAPH_COMBINED_VCFS) \
	nanopore_sv_calling/SV_NORMALIZATION $(NANOPORE_NORMALIZED_SVS) \
	sv_genotyping/combined_svs/sveval_benchmarks/benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/combined_svs/sveval_benchmarks ; $(R_RUN_COMMAND) benchmark.R

# Benchmark of combined Illumina/Oxford Nanopore SVs in non-repeat regions
sv_genotyping/combined_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData : \
	sv_genotyping/combined_svs/PARAGRAPH_COMBINED_GENOTYPING $(PARAGRAPH_COMBINED_VCFS) \
	nanopore_sv_calling/SV_NORMALIZATION $(NANOPORE_NORMALIZED_SVS) \
	sv_genotyping/combined_svs/sveval_benchmarks/norepeat_benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/repeat_regions/non_repeated_regions.bed \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/combined_svs/sveval_benchmarks ; $(R_RUN_COMMAND) norepeat_benchmark.R

# --- This section generates the results for the analysis of refined versus raw SVs

# Filtering and normalizing the raw SVs used for the benchmarks
breakpoint_refinement_analysis/raw_svs/RAW_SV_CALLS : \
	nanopore_sv_calling/SV_CALLING $(NANOPORE_FILTERED_SVS) \
	breakpoint_refinement_analysis/raw_svs/process_vcf_files.sh \
	breakpoint_refinement_analysis/raw_svs/fix_vcfs.R \
	breakpoint_refinement_analysis/nanopore_ids.txt \
	scripts/fix_sniffles.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd breakpoint_refinement_analysis/raw_svs ; ./process_vcf_files.sh $(BCFTOOLS) $(R_RUN_COMMAND) ; touch RAW_SV_CALLS

# Merging the raw variants together with SVmerge
breakpoint_refinement_analysis/raw_svs/svmerged.clustered.vcf : breakpoint_refinement_analysis/raw_svs/RAW_SV_CALLS \
	breakpoint_refinement_analysis/raw_svs/SVmerge_variants.sh \
	breakpoint_refinement_analysis/raw_svs/files.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd breakpoint_refinement_analysis/raw_svs ; ./SVmerge_variants.sh $(SVMERGE)

# Merging the refined variants together with SVmerge
breakpoint_refinement_analysis/refined_svs/svmerged_clustered_sorted.vcf : nanopore_sv_calling/SV_NORMALIZATION \
	breakpoint_refinement_analysis/refined_svs/SVmerge_variants.sh \
	breakpoint_refinement_analysis/refined_svs/files.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	breakpoint_refinement_analysis/refined_svs/select_svs.R \
	scripts/merge_realigned_variants.R \
	breakpoint_refinement_analysis/refined_svs/process_sort_svmerged.sh
	cd breakpoint_refinement_analysis/refined_svs ; ./SVmerge_variants.sh $(SVMERGE) ; process_sort_svmerged.sh $(R_RUN_COMMAND) $(BCFTOOLS)

# Computing the k-mers for BayesTyper
breakpoint_refinement_analysis/bayestyper_kmers/BAYESTYPER_KMERS : illumina_data/ILLUMINA_TRIMMING \
	breakpoint_refinement_analysis/bayestyper_kmers/kmc_bloom.sh \
	breakpoint_refinement_analysis/illumina_ids.txt
	cd breakpoint_refinement_analysis/bayestyper_kmers ; ./kmc_bloom.sh $(KMC) $(BAYESTYPERTOOLS); touch BAYESTYPER_KMERS

# Genotyping raw variants with Bayestyper
breakpoint_refinement_analysis/raw_svs/bayestyper/bayestyper_split.vcf : \
	breakpoint_refinement_analysis/raw_svs/svmerged.clustered.vcf \
	breakpoint_refinement_analysis/bayestyper_kmers/BAYESTYPER_KMERS \
	breakpoint_refinement_analysis/raw_svs/bayestyper/bayestyper.sh \
	breakpoint_refinement_analysis/raw_svs/bayestyper/samples.tsv \
	refgenome/Gmax_508_v4.0.fa \
	refgenome/bt_decoy_sequences.fasta
	cd breakpoint_refinement_analysis/raw_svs/bayestyper ; ./bayestyper.sh $(BAYESTYPER) $(BAYESTYPERTOOLS) $(BCFTOOLS)

# Genotyping raw variants with Paragraph
breakpoint_refinement_analysis/raw_svs/paragraph/PARAGRAPH_RAW_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT $(ILLUMINA_ALIGNED_READS) \
	sv_genotyping/MANIFEST_FILES $(PARAGRAPH_MANIFEST_FILES) \
	breakpoint_refinement_analysis/raw_svs/svmerged.clustered.vcf \
	breakpoint_refinement_analysis/raw_svs/paragraph/run_paragraph.sh \
	scripts/addMissingPaddingGmax4.py \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	breakpoint_refinement_analysis/illumina_ids.txt
	cd breakpoint_refinement_analysis/raw_svs/paragraph ; ./run_paragraph.sh $(BCFTOOLS) $(BGZIP) $(TABIX) $(MULTIGRMPY) ; touch PARAGRAPH_RAW_GENOTYPING

# Genotyping raw variants with vg
breakpoint_refinement_analysis/raw_svs/vg/VG_RAW_GENOTYPING : illumina_data/ILLUMINA_TRIMMING \
	breakpoint_refinement_analysis/raw_svs/vg/vg.sh \
	breakpoint_refinement_analysis/raw_svs/svmerged.clustered.vcf \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	breakpoint_refinement_analysis/illumina_ids.txt 
	cd breakpoint_refinement_analysis/raw_svs/vg ; ./vg.sh $(BGZIP) $(TABIX) $(VG) ; touch VG_RAW_GENOTYPING

# Benchmarking the raw variants genotyped with all three genotyping tools
breakpoint_refinement_analysis/raw_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData : \
	breakpoint_refinement_analysis/raw_svs/RAW_SV_CALLS \
	breakpoint_refinement_analysis/raw_svs/sveval_benchmarks/raw_nogeno_analysis.R \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	breakpoint_refinement_analysis/raw_svs/bayestyper/bayestyper_split.vcf \
	breakpoint_refinement_analysis/raw_svs/paragraph/PARAGRAPH_RAW_GENOTYPING \
	breakpoint_refinement_analysis/raw_svs/vg/VG_RAW_GENOTYPING
	cd breakpoint_refinement_analysis/raw_svs/sveval_benchmarks ; $(R_RUN_COMMAND) raw_nogeno_analysis.R

# Genotyping refined variants with Bayestyper
breakpoint_refinement_analysis/refined_svs/bayestyper/bayestyper_split.vcf : \
	breakpoint_refinement_analysis/refined_svs/svmerged_clustered_sorted.vcf \
	breakpoint_refinement_analysis/bayestyper_kmers/BAYESTYPER_KMERS \
	breakpoint_refinement_analysis/refined_svs/bayestyper/bayestyper.sh \
	breakpoint_refinement_analysis/refined_svs/bayestyper/samples.tsv \
	refgenome/Gmax_508_v4.0.fa \
	refgenome/bt_decoy_sequences.fasta
	cd breakpoint_refinement_analysis/refined_svs/bayestyper ; ./bayestyper.sh $(BAYESTYPER) $(BAYESTYPERTOOLS) $(BCFTOOLS)

# Genotyping refined variants with Paragraph
breakpoint_refinement_analysis/refined_svs/paragraph/PARAGRAPH_REFINED_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT $(ILLUMINA_ALIGNED_READS) \
	sv_genotyping/MANIFEST_FILES $(PARAGRAPH_MANIFEST_FILES) \
	breakpoint_refinement_analysis/refined_svs/svmerged_clustered_sorted.vcf \
	breakpoint_refinement_analysis/refined_svs/paragraph/run_paragraph.sh \
	scripts/addMissingPaddingGmax4.py \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	breakpoint_refinement_analysis/illumina_ids.txt
	cd breakpoint_refinement_analysis/refined_svs/paragraph ; ./run_paragraph.sh $(BCFTOOLS) $(BGZIP) $(TABIX) $(MULTIGRMPY) ; touch PARAGRAPH_REFINED_GENOTYPING

# Genotyping refined variants with vg
breakpoint_refinement_analysis/refined_svs/vg/VG_REFINED_GENOTYPING : illumina_data/ILLUMINA_TRIMMING \
	breakpoint_refinement_analysis/refined_svs/vg/vg.sh \
	breakpoint_refinement_analysis/refined_svs/svmerged_clustered_sorted.vcf \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	breakpoint_refinement_analysis/illumina_ids.txt 
	cd breakpoint_refinement_analysis/refined_svs/vg ; ./vg.sh $(BGZIP) $(TABIX) $(VG) ; touch VG_REFINED_GENOTYPING

# Benchmarking the refined variants genotyped with all three genotyping tools
breakpoint_refinement_analysis/refined_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData : \
	nanopore_sv_calling/SV_NORMALIZATION \
	breakpoint_refinement_analysis/refined_svs/sveval_benchmarks/refined_nogeno_analysis.R \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	breakpoint_refinement_analysis/refined_svs/bayestyper/bayestyper_split.vcf \
	breakpoint_refinement_analysis/refined_svs/paragraph/PARAGRAPH_REFINED_GENOTYPING \
	breakpoint_refinement_analysis/refined_svs/vg/VG_REFINED_GENOTYPING
	cd breakpoint_refinement_analysis/refined_svs/sveval_benchmarks ; $(R_RUN_COMMAND) refined_nogeno_analysis.R

# --- This section performs the fastStructure analysis on SNPs

# Calling the SNPs using Platypus
structure_analysis/platypus_snps.vcf : illumina_data/ILLUMINA_ALIGNMENT \
	structure_analysis/call_snps.sh \
	structure_analysis/platypus102_Gmax_v4_bamfiles.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd structure_analysis ; ./call_snps.sh $(PLATYPUS)

# Filtering the SNPs called with Platypus
structure_analysis/platypus_filtered_snps.vcf : structure_analysis/platypus_snps.vcf \
	structure_analysis/filter_platypus_snps.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/het_counts.awk \
	scripts/het_stats.awk \
	scripts/filter_hetsites.R
	cd structure_analysis ; ./filter_platypus_snps.sh $(BCFTOOLS) $(VCFTOOLS) $(R_RUN_COMMAND)

# Computing the fastStructure results from the filtered Platypus SNPs
structure_analysis/structure.5.meanQ : structure_analysis/platypus_filtered_snps.vcf \
	structure_analysis/compute_structure.sh
	cd structure_analysis ; ./compute_structure.sh $(FASTSTRUCTURE) $(VCFTOOLS) $(PLINK)

# --- This section analyses the TEs found in the combined Illumina/Oxford Nanopore SV dataset
te_analysis/polymorphic_tes.tsv : te_analysis/query_all.vcf \
	te_analysis/blast_svs.txt \
	te_analysis/te_blast_analysis.R \
	te_analysis/extract_te_metadata.sh \
	te_analysis/te_database/SoyBase_TE_Fasta.txt
	cd te_analysis ; $(R_RUN_COMMAND) te_blast_analysis.R

te_analysis/query_all.vcf : sv_genotyping/combined_svs/combined_paragraph_filtered.vcf \
	te_analysis/extract_query_vcf.sh
	cd te_analysis ; ./extract_query_vcf.sh

te_analysis/blast_svs.txt : te_analysis/query_all.vcf \
	te_analysis/blast_tes.sh \
	te_analysis/te_database/SoyBase_TE_Fasta.txt
	cd te_analysis ; ./blast_tes.sh $(BLASTN) $(MAKEBLASTDB)

