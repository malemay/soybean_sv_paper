
##### CREATING SOME VARIABLES FOR EXECUTABLES

AGE = ~/programs/AGE/age_align
ASMVAR = ~/programs/AsmVar/src/AsmvarDetect/ASV_VariantDetector
BAYESTYPER = ~/programs/bayesTyper_v1.5_linux_x86_64/bin/bayesTyper
BAYESTYPERTOOLS = ~/programs/bayesTyper_v1.5_linux_x86_64/bin/bayesTyperTools
BBDUK = ~/programs/bbmap/bbduk.sh
BAMADDRG = /prg/bamaddrg/1.0/bamaddrg
BGZIP = /prg/htslib/1.10.2/bin/bgzip
BWA = /prg/bwa/0.7.17/bwa
BCFTOOLS = ~/programs/bcftools/bcftools
BCFTOOLS_PLUGIN_PATH = ~/programs/bcftools/plugins
BLASTN = ~/programs/ncbi-blast-2.11.0+/bin/blastn
CIRCOS = ~/programs/circos-0.69-9/bin/circos
FASTSTRUCTURE = /prg/fastStructure/1.0/structure.py
FLASH = ~/programs/FLASH-1.2.11-Linux-x86_64/flash
GINSI = ~/.local/bin/ginsi
GRF = ~/programs/GenericRepeatFinder/bin/grf-main
IDXDEPTH = ~/programs/paragraph/bin/idxdepth
KMC = ~/programs/KMC3.linux/kmc
LASTAL = ~/programs/last-1047/src/lastal
LASTSPLIT = ~/programs/last-1047/src/last-split
MAKEBLASTDB = ~/programs/ncbi-blast-2.11.0+/bin/makeblastdb
MANTA = ~/programs/manta-1.6.0.centos6_x86_64/bin/configManta.py
MINIMAP2 = ~/programs/minimap2/minimap2
MULTIGRMPY = ~/programs/paragraph/bin/multigrmpy.py
NANOPLOT = ~/.local/bin/NanoPlot
NGMLR = ~/programs/ngmlr/ngmlr-0.2.7/ngmlr
PLATYPUS = /prg/platypus/0.8.1.1/bin/Platypus.py
PLINK = /prg/plink/1.90b5.3/bin/plink
PORECHOP = ~/programs/Porechop/porechop-runner.py
# Using R 3.5.0 for figures because of a problem with resolution when using R 4.0
R_FIG_COMMAND = /usr/bin/Rscript
R_RUN_COMMAND = /prg/R/4.0/bin/Rscript
SAMTOOLS = ~/programs/samtools/samtools
SMOOVE = ~/programs/smoove/smoove
SOAPDENOVO2 = /prg/SOAPdenovo/2.04/SOAPdenovo-63mer
SNIFFLES = ~/programs/Sniffles-master/bin/sniffles-core-1.0.11/sniffles
SVABA = ~/programs/svaba/bin/svaba
SVMERGE = ~/programs/SVanalyzer-install/bin/SVmerge
TABIX = /prg/htslib/1.10.2/bin/tabix
VCFTOOLS = /prg/vcftools/0.1.16/bin/vcftools
VG = ~/programs/vg/vg
WTDBG2 = ~/programs/wtdbg2/wtdbg2
WTPOA_CNS = ~/programs/wtdbg2/wtpoa-cns


##### CREATING SOME VARIABLES AND MAIN TARGETS

# Figures 1 to 6
FIGURES := $(shell seq 1 6 | xargs -I {} echo figures/figure_{}.png)

# Tables 1 to 3
TABLES := $(shell seq 1 3 | xargs -I {} echo tables/table_{}.png)

# Supplemental files
SDIR := additional_files
SUPFILES := $(SDIR)/additional_file_1.pdf $(SDIR)/additional_te_file.csv \
	$(SDIR)/additional_bp_over_file.csv $(SDIR)/additional_bp_under_file.csv \
	$(SDIR)/additional_pfam_over_file.csv  $(SDIR)/additional_pfam_under_file.csv \
	$(SDIR)/additional_tir_similarity_file.csv

# Supplemental figures S1 to S19 plus figure S20 which depends on another file
SUPFIGURES := $(shell seq 1 19 | xargs -I {} echo figures/figure_s{}.png) figures/Gm04_2257090_annotated.png 

# Supplemental tables S1 to S8
SUPTABLES := $(shell seq 1 8 | xargs -I {} echo tables/table_s{}.csv)

# --- This target prepares all the figures, tables, and supplemental data
all: $(FIGURES) $(TABLES) $(SUPFILES)

# A target for all main figures
figures : $(FIGURES)

# A target for all main tables
tables : $(TABLES)

# A target for all supplemental figures
supfigures : $(SUPFIGURES)

# A target for all supplemental tables
subtables : $(SUPTABLES)

# A target for all supplemental files
supfiles : $(SUPFILES)

# Compiling the Supplemental Data file from the .tex file as well as supplementary tables and figures
$(SDIR)/additional_file_1.pdf : $(SUPFIGURES) $(SUPTABLES) \
	$(SDIR)/additional_file_1.tex $(SDIR)/references.bib $(SDIR)/vancouver.bst
	cd $(SDIR) ; pdflatex additional_file_1.tex; bibtex additional_file_1; pdflatex additional_file_1.tex; pdflatex additional_file_1.tex

# Generating the additional CSV files
$(SDIR)/additional_%.csv : $(SDIR)/additional_%.R
	cd $(SDIR) ; $(R_RUN_COMMAND) $(<F)

$(SDIR)/additional_te_file.csv : te_analysis/polymorphic_tes.tsv
$(SDIR)/additional_bp_over_file.csv : scripts/format_go_csv.R gene_analysis/GO_ANALYSIS
$(SDIR)/additional_bp_under_file.csv : scripts/format_go_csv.R gene_analysis/GO_ANALYSIS
$(SDIR)/additional_pfam_over_file.csv : scripts/format_go_csv.R gene_analysis/GO_ANALYSIS
$(SDIR)/additional_pfam_under_file.csv : scripts/format_go_csv.R gene_analysis/GO_ANALYSIS
$(SDIR)/additional_tir_similarity_file.csv : te_analysis/multiple_alignments/TIR_TSD_ANALYSIS

##### CREATING THE MAIN FIGURES

# Creating a pattern rule that will be used to run the R commands for every figure but figure 3
figures/figure_%.png : figures/figure_%.R
	cd figures; $(R_FIG_COMMAND) $(<F)

figures/figure_1.png : sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData scripts/make_plot_data.R
figures/figure_2.png : sv_genotyping/nanopore_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData scripts/make_plot_data.R 

# --- Section for the Circos figure
CIRCD := figures/figure_3_circos
EXTDIR := external/circos_config_files

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

figures/figure_4.png : structure_analysis/snp_pca/SNP_PCA structure_analysis/sv_pca/SV_PCA structure_analysis/structure.5.meanQ
figures/figure_5.png : gene_analysis/GENE_OVERLAP_ANALYSIS gene_analysis/permutation_all_100kb.RData

figures/figure_6.png : te_analysis/polymorphic_tes.tsv \
	te_analysis/tian2012_tes.txt \
	te_analysis/multiple_alignments/TIR_TSD_ANALYSIS \
	te_analysis/multiple_alignments/Gm04_2257090_INS_480_analysis/STOWAWAY_MITE_ANALYSIS


#### CREATING THE MAIN RESULTS TABLES

# Creating a pattern rule that will be used to generate the tables in PNG format from the R code and .tex file
tables/%.png : tables/%.R tables/%.tex
	cd tables ; $(R_RUN_COMMAND) $(*F).R ; pdflatex $(*F).tex ; pdftoppm $(*F).pdf -singlefile -r 500 -png $(*F)

tables/table_1.png : tables/gather_table_1_data.sh \
	illumina_sv_calling/asmvar/asmvar_filtering/asmvar_svs.vcf \
	illumina_sv_calling/manta/manta_svs.vcf \
	illumina_sv_calling/smoove/smoove_svs.vcf \
	illumina_sv_calling/svaba/svaba_svs.vcf \
	sv_genotyping/illumina_svs/svmerged.clustered.vcf \
	scripts/count_svtypes_svsizes.awk

tables/table_2.png : te_analysis/polymorphic_tes.tsv refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3
tables/table_3.png : gene_analysis/GENE_OVERLAP_ANALYSIS


##### CREATING THE SUPPLEMENTAL FIGURES

# These use the same pattern rule as the main figures so we do not need to specify it again

figures/figure_s1.png : sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData scripts/make_plot_data.R
figures/figure_s2.png : sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData scripts/make_plot_data.R
figures/figure_s3.png : sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData scripts/make_plot_data.R

figures/figure_s4.png : depth_distributions/average_depth.RData utilities/line_ids.txt \
	sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData

figures/figure_s5.png : sv_genotyping/illumina_svs/sveval_benchmarks/frequency_RData/sveval_frequency_rates.RData scripts/make_plot_data.R
figures/figure_s6.png : sv_genotyping/illumina_svs/sveval_benchmarks/frequency_RData/sveval_frequency_rates.RData scripts/make_plot_data.R
figures/figure_s7.png : sv_genotyping/illumina_svs/sveval_benchmarks/NCALLERS_ILLUMINA_BENCHMARK scripts/make_plot_data.R
figures/figure_s8.png : sv_genotyping/illumina_svs/size_distribution.tsv
figures/figure_s9.png : sv_genotyping/illumina_svs/sveval_benchmarks/NCALLERS_ILLUMINA_BENCHMARK scripts/format_sveval_plotting_data.R
figures/figure_s10.png : sv_genotyping/illumina_svs/sveval_benchmarks/NCALLERS_ILLUMINA_BENCHMARK scripts/format_sveval_plotting_data.R

figures/figure_s11.png : nanoplot_stats/N50_stats.txt utilities/line_ids.txt \
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


##### CREATING THE SUPPLEMENTAL TABLES

# Creating a pattern rule that generates the tables in csv format from the R code
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


##### CREATING SOME INTERMEDIATE DATA FOR THE FIGURES AND TABLES

# Computing the average sequencing depths used for plotting figure S4
depth_distributions/average_depth.RData : nanopore_data/NANOPORE_ALIGNMENT \
	utilities/line_ids.txt \
	depth_distributions/depth_distributions.sh \
	scripts/depth_distribution.sh \
	depth_distributions/compute_average_depth.R
	cd depth_distributions/ ; ./depth_distributions.sh $(SAMTOOLS) ; $(R_RUN_COMMAND) compute_average_depth.R

# Generating the SV size distributions necessary for figure S8
sv_genotyping/illumina_svs/size_distribution.tsv : sv_genotyping/illumina_svs/svmerged.clustered.vcf \
	sv_genotyping/illumina_svs/extract_size_distribution.sh
	cd sv_genotyping/illumina_svs/ ; ./extract_size_distribution.sh $(BCFTOOLS)

# Generating the Nanoplot N50 data for figure S11
nanoplot_stats/NANOPLOT_BAM_STATS : nanopore_data/NANOPORE_ALIGNMENT \
	nanoplot_stats/nanoplot_sorted_bam_all.sh \
	utilities/line_ids.txt
	cd nanoplot_stats ; ./nanoplot_sorted_bam_all.sh $(NANOPLOT) ; touch nanoplot_stats/NANOPLOT_BAM_STATS

nanoplot_stats/N50_stats.txt : nanoplot_stats/NANOPLOT_BAM_STATS \
	nanoplot_stats/gather_N50_stats.sh \
	utilities/line_ids.txt
	cd nanoplot_stats ; ./gather_N50_stats.sh > N50_stats.txt

# Generating the data on the number of SVs per sample for Table S3
tables/table_s3_data.txt : nanopore_sv_calling/SV_NORMALIZATION \
       	tables/gather_table_s3_data.sh \
       	utilities/line_ids.txt \
	scripts/count_svtypes_svsizes.awk
	cd tables ; ./gather_table_s3_data.sh

# Generating the data on the testing of the breakpoint refinement pipeline for Table S4
nanopore_sv_calling/all_metainfo.RData : nanopore_sv_calling/SV_NORMALIZATION \
	nanopore_sv_calling/gather_metainfo.R
	cd nanopore_sv_calling/ ; $(R_RUN_COMMAND) gather_metainfo.R $(BCFTOOLS)


##### ANALYZING THE OXFORD NANOPORE SEQUENCING DATA

# This variable stores the paths to the Oxford Nanopore fastq files
NANOPORE_READS := $(shell cat utilities/flowcell_names.txt | xargs -I {} echo nanopore_data/{}.fastq.gz)

# Processing the raw basecalled Nanopore reads using Porechop and aligning them to the reference genome
nanopore_data/NANOPORE_ALIGNMENT : $(NANOPORE_READS) \
	nanopore_data/read_processing_alignment.sh \
	utilities/flowcell_names.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd nanopore_data ; ./read_processing_alignment.sh $(PORECHOP) $(NGMLR) $(SAMTOOLS) ; touch NANOPORE_ALIGNMENT

# Calling SVs with Sniffles and filtering the output
nanopore_sv_calling/SV_CALLING : nanopore_data/NANOPORE_ALIGNMENT \
	nanopore_sv_calling/sniffles_calling.sh \
	nanopore_sv_calling/filter_sniffles_vcfs.R \
	scripts/filter_sniffles.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd nanopore_sv_calling ; ./sniffles_calling.sh $(SNIFFLES) $(R_RUN_COMMAND) ; touch SV_CALLING

# Refining the SV breakpoints
nanopore_sv_calling/SV_REFINEMENT : nanopore_sv_calling/SV_CALLING \
	nanopore_data/NANOPORE_ALIGNMENT \
	scripts/breakpoint_refinement/breakpoint_refinement.R \
	utilities/line_ids.txt \
	scripts/breakpoint_refinement/age_realign.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd nanopore_sv_calling ; $(R_RUN_COMMAND) refine_breakpoints.R $(SAMTOOLS) $(MINIMAP2) $(AGE) $(WTDBG2) $(WTPOA_CNS) ; touch SV_REFINEMENT

# Processing the refined VCFs to prepare them for input to SVmerge;
# These files are also the reference Oxford Nanpore SVs for benchmarking
nanopore_sv_calling/SV_NORMALIZATION : nanopore_sv_calling/SV_REFINEMENT \
	nanopore_sv_calling/process_vcf_files.sh \
	nanopore_sv_calling/add_metainfo_all.R \
	nanopore_sv_calling/fix_vcfs.R \
	scripts/fix_sniffles.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd nanopore_sv_calling ; ./process_vcf_files.sh $(BCFTOOLS) $(R_RUN_COMMAND) ; touch SV_NORMALIZATION


##### ANALYZING THE ILLUMINA SEQUENCING DATA

# These variables store the paths to the location of the Illumina reads
ILLUMINA_READS_1 := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_data/raw_fastq/{}_R1.fastq.gz)
ILLUMINA_READS_2 := $(shell cat utilities/all_lines.txt | xargs -I {} echo illumina_data/raw_fastq/{}_R2.fastq.gz)

# Trimming the Illumina reads with bbduk
illumina_data/ILLUMINA_TRIMMING : $(ILLUMINA_READS_1) $(ILLUMINA_READS_2) \
	illumina_data/bbduk_trimming.sh \
	utilities/all_lines.txt
	cd illumina_data ; ./bbduk_trimming.sh $(BBDUK) ; touch ILLUMINA_TRIMMING

# Aligning the trimmed Illumina reads with bwa mem
illumina_data/ILLUMINA_ALIGNMENT : illumina_data/ILLUMINA_TRIMMING \
	illumina_data/bwa_alignment.sh \
	utilities/all_lines.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd illumina_data ; ./bwa_alignment.sh $(BWA) $(SAMTOOLS) $(BAMADDRG) ; touch ILLUMINA_ALIGNMENT


##### CALLING THE SVS ON THE ILLUMINA DATA

# Merging the trimmed sequences with FLASH
illumina_sv_calling/asmvar/FLASH_MERGING : illumina_data/ILLUMINA_TRIMMING  \
	illumina_sv_calling/asmvar/flash_merging.sh
	cd illumina_sv_calling/asmvar ; ./flash_merging.sh $(FLASH) ; touch FLASH_MERGING

# Assembling the sequences using SOAPdenovo2
illumina_sv_calling/asmvar/SOAPDENOVO_ASSEMBLY : illumina_data/ILLUMINA_TRIMMING \
	illumina_sv_calling/asmvar/FLASH_MERGING \
	illumina_sv_calling/asmvar/assembly.sh \
	illumina_sv_calling/asmvar/soap_assembly/make_config.sh \
	illumina_sv_calling/asmvar/soap_assembly/config_template.txt \
	utilities/all_lines.txt
	cd illumina_sv_calling/asmvar ; ./assembly.sh $(SOAPDENOVO2) ; touch SOAPDENOVO_ASSEMBLY

# Aligning the de novo assemblies to the reference genome using LAST
illumina_sv_calling/asmvar/LAST_ALIGNMENT : illumina_sv_calling/asmvar/SOAPDENOVO_ASSEMBLY \
	illumina_sv_calling/asmvar/last_alignment.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd illumina_sv_calling/asmvar ; ./last_alignment.sh $(LASTAL) $(LASTSPLIT) ; touch LAST_ALIGNMENT

# Calling SVs with AsmVar on the aligned assemblies
illumina_sv_calling/asmvar/ASMVAR_CALLING : illumina_sv_calling/asmvar/LAST_ALIGNMENT \
	illumina_sv_calling/asmvar/SOAPDENOVO_ASSEMBLY \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd illumina_sv_calling/asmvar ; ./asmvar_call.sh $(ASMVAR) ; touch ASMVAR_CALLING

# Filtering the VCFs resulting from calling AsmVar on assemblies
illumina_sv_calling/asmvar/asmvar_filtering/asmvar_svs.vcf : illumina_sv_calling/asmvar/ASMVAR_CALLING \
	illumina_sv_calling/asmvar/asmvar_filter.sh \
	utilities/all_lines.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai \
	scripts/add_svtype.awk \
	illumina_sv_calling/asmvar/svtype_header_line.txt \
	scripts/extract_svs_50.awk
	cd illumina_sv_calling/asmvar ; ./asmvar_filter.sh $(BCFTOOLS) $(BAYESTYPERTOOLS)

# Calling the SVs using manta
illumina_sv_calling/manta/MANTA_CALLING : illumina_data/ILLUMINA_ALIGNMENT \
	illumina_sv_calling/manta/manta_call.sh \
	illumina_sv_calling/manta/manta_config.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd illumina_sv_calling/manta ; ./manta_call.sh $(MANTA) ; touch MANTA_CALLING

# Filtering the SVs discovered using Manta
illumina_sv_calling/manta/manta_svs.vcf : illumina_sv_calling/manta/MANTA_CALLING \
	illumina_sv_calling/manta/manta_filter.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	illumina_sv_calling/manta/ACO_header_line.txt \
	scripts/extract_svs_50.awk
	cd illumina_sv_calling/manta ; ./manta_filter.sh $(BCFTOOLS) $(BAYESTYPERTOOLS) $(TABIX)

# Calling the SVs using SvABA
illumina_sv_calling/svaba/SVABA_CALLING : illumina_data/ILLUMINA_ALIGNMENT \
	illumina_sv_calling/svaba/svaba_call.sh \
	utilities/all_lines.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd illumina_sv_calling/svaba ; ./svaba_call.sh $(SVABA) ; touch SVABA_CALLING

# Filtering the SVs discovered using SvABA
illumina_sv_calling/svaba/svaba_svs.vcf : illumina_sv_calling/svaba/SVABA_CALLING \
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

# Calling the SVs using smoove
illumina_sv_calling/smoove/all_samples.smoove.square.vcf.gz : illumina_data/ILLUMINA_ALIGNMENT \
	illumina_sv_calling/smoove/smoove_call.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd illumina_sv_calling/smoove ; ./smoove_call.sh $(SMOOVE) $(TABIX)

# Filtering the SVs called using smoove
illumina_sv_calling/smoove/smoove_svs.vcf : illumina_sv_calling/smoove/all_samples.smoove.square.vcf.gz \
	illumina_sv_calling/smoove/smoove_filter.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	illumina_sv_calling/smoove/ACO_header_line.txt \
	scripts/extract_svs_50.awk
	cd illumina_sv_calling/smoove ; ./smoove_filter.sh $(BCFTOOLS) $(BAYESTYPERTOOLS)


##### GENOTYPING THE VARIOUS SV DATASETS USING PARAGRAPH

# Preparing the manifest files needed by Paragraph, irrespective of the dataset to genotype
sv_genotyping/MANIFEST_FILES : illumina_data/ILLUMINA_ALIGNMENT \
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

sv_genotyping/illumina_svs/PARAGRAPH_ILLUMINA_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT \
	sv_genotyping/MANIFEST_FILES \
	sv_genotyping/illumina_svs/svmerged.clustered.vcf \
	sv_genotyping/illumina_svs/run_paragraph.sh \
	scripts/addMissingPaddingGmax4.py \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/all_lines.txt
	cd sv_genotyping/illumina_svs ; ./run_paragraph.sh $(BCFTOOLS) $(BGZIP) $(TABIX) $(MULTIGRMPY) ; touch PARAGRAPH_ILLUMINA_GENOTYPING

# --- This section genotypes the Oxford Nanopore SVs using Paragraph, first preparing them with SVmerge
sv_genotyping/nanopore_svs/svmerged_clustered_sorted.vcf : sv_genotyping/nanopore_svs/svmerge_files.txt \
	nanopore_sv_calling/SV_NORMALIZATION \
	sv_genotyping/nanopore_svs/SVmerge.sh \
	sv_genotyping/nanopore_svs/merge_realigned.R \
	scripts/merge_realigned_variants.R \
	sv_genotyping/nanopore_svs/svmerged_sort.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd sv_genotyping/nanopore_svs ; ./SVmerge.sh $(SVMERGE) ; $(R_RUN_COMMAND) merge_realigned.R ; ./svmerged_sort.sh $(BCFTOOLS)

sv_genotyping/nanopore_svs/PARAGRAPH_NANOPORE_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT \
	sv_genotyping/MANIFEST_FILES \
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

sv_genotyping/combined_svs/PARAGRAPH_COMBINED_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT \
	sv_genotyping/MANIFEST_FILES \
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


##### BENCHMARKING THE SV CALLS

# Benchmark of Illumina SVs genome-wide
sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData : \
	sv_genotyping/illumina_svs/PARAGRAPH_ILLUMINA_GENOTYPING \
	nanopore_sv_calling/SV_NORMALIZATION \
	sv_genotyping/illumina_svs/sveval_benchmarks/benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; $(R_RUN_COMMAND) benchmark.R

# Benchmark of Illumina SVs in non-repeat regions
sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData : \
	sv_genotyping/illumina_svs/PARAGRAPH_ILLUMINA_GENOTYPING \
	nanopore_sv_calling/SV_NORMALIZATION \
	sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/repeat_regions/non_repeated_regions.bed \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; $(R_RUN_COMMAND) norepeat_benchmark.R

# Benchmark of the Illumina SVs as a function of the homozygous ALT count
# For this we first need to prepare the input vcf file for the benchmark
sv_genotyping/illumina_svs/illumina_paragraph_merged.vcf : sv_genotyping/illumina_svs/PARAGRAPH_ILLUMINA_GENOTYPING \
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
	nanopore_sv_calling/SV_NORMALIZATION \
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
	nanopore_sv_calling/SV_NORMALIZATION \
	sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; $(R_RUN_COMMAND) ncallers_benchmark.R

# Benchmark of Oxford Nanopore SVs genome-wide
sv_genotyping/nanopore_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData : \
	sv_genotyping/nanopore_svs/PARAGRAPH_NANOPORE_GENOTYPING \
	nanopore_sv_calling/SV_NORMALIZATION \
	sv_genotyping/nanopore_svs/sveval_benchmarks/benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/nanopore_svs/sveval_benchmarks ; $(R_RUN_COMMAND) benchmark.R

# Benchmark of Oxford Nanopore SVs in non-repeat regions
sv_genotyping/nanopore_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData : \
	sv_genotyping/nanopore_svs/PARAGRAPH_NANOPORE_GENOTYPING \
	nanopore_sv_calling/SV_NORMALIZATION \
	sv_genotyping/nanopore_svs/sveval_benchmarks/norepeat_benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/repeat_regions/non_repeated_regions.bed \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/nanopore_svs/sveval_benchmarks ; $(R_RUN_COMMAND) norepeat_benchmark.R

# Benchmark of combined Illumina/Oxford Nanopore SVs genome-wide
sv_genotyping/combined_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData : \
	sv_genotyping/combined_svs/PARAGRAPH_COMBINED_GENOTYPING \
	nanopore_sv_calling/SV_NORMALIZATION \
	sv_genotyping/combined_svs/sveval_benchmarks/benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/combined_svs/sveval_benchmarks ; $(R_RUN_COMMAND) benchmark.R

# Benchmark of combined Illumina/Oxford Nanopore SVs in non-repeat regions
sv_genotyping/combined_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData : \
	sv_genotyping/combined_svs/PARAGRAPH_COMBINED_GENOTYPING \
	nanopore_sv_calling/SV_NORMALIZATION \
	sv_genotyping/combined_svs/sveval_benchmarks/norepeat_benchmark.R \
	utilities/line_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/repeat_regions/non_repeated_regions.bed \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/combined_svs/sveval_benchmarks ; $(R_RUN_COMMAND) norepeat_benchmark.R


##### ANALYSIS OF THE RAW VERSUS REFINED SVS (FOR SUPPLEMENTARY MATERIAL)

# Filtering and normalizing the raw SVs used for the benchmarks
breakpoint_refinement_analysis/raw_svs/RAW_SV_CALLS : \
	nanopore_sv_calling/SV_CALLING \
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
breakpoint_refinement_analysis/raw_svs/paragraph/PARAGRAPH_RAW_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT \
	sv_genotyping/MANIFEST_FILES \
	breakpoint_refinement_analysis/raw_svs/svmerged.clustered.vcf \
	breakpoint_refinement_analysis/raw_svs/paragraph/run_paragraph.sh \
	scripts/addMissingPaddingGmax4.py \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	breakpoint_refinement_analysis/illumina_ids.txt
	cd breakpoint_refinement_analysis/raw_svs/paragraph ; ./run_paragraph.sh $(BCFTOOLS) $(BGZIP) $(TABIX) $(MULTIGRMPY) ; \
		touch PARAGRAPH_RAW_GENOTYPING

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
breakpoint_refinement_analysis/refined_svs/paragraph/PARAGRAPH_REFINED_GENOTYPING : illumina_data/ILLUMINA_ALIGNMENT \
	sv_genotyping/MANIFEST_FILES \
	breakpoint_refinement_analysis/refined_svs/svmerged_clustered_sorted.vcf \
	breakpoint_refinement_analysis/refined_svs/paragraph/run_paragraph.sh \
	scripts/addMissingPaddingGmax4.py \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	breakpoint_refinement_analysis/illumina_ids.txt
	cd breakpoint_refinement_analysis/refined_svs/paragraph ; ./run_paragraph.sh $(BCFTOOLS) $(BGZIP) $(TABIX) $(MULTIGRMPY) ; \
		touch PARAGRAPH_REFINED_GENOTYPING

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


##### POPULATION STRUCTURE ANALYSES FOR FIGURE 4

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

# Computing the PCA on the SNPs called with Platypus
structure_analysis/snp_pca/SNP_PCA : structure_analysis/platypus_filtered_snps.vcf \
	structure_analysis/snp_pca/snp_pca.sh
	cd structure_analysis/snp_pca ; ./snp_pca.sh $(VCFTOOLS) $(PLINK) ; touch SNP_PCA

# Computing the PCA on the combined Illumina/Oxford Nanopore SVs genotyped by Paragraph
structure_analysis/sv_pca/SV_PCA : sv_genotyping/combined_svs/combined_paragraph_filtered.vcf \
	structure_analysis/sv_pca/sv_pca.sh \
	scripts/recode_alleles.awk
	cd structure_analysis/sv_pca ; ./sv_pca.sh $(VCFTOOLS) $(PLINK) ; touch SV_PCA


##### GENE FEATURE OVERLAP ANALYSES

# Creating the bed file of repeats from the reference genome and Phytozome repeat annotation
refgenome/repeat_regions/non_repeated_regions.bed : refgenome/repeat_regions/make_repeat_bed.R \
	refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3 \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd refgenome/repeat_regions ; $(R_RUN_COMMAND) make_repeat_bed.R

# Generating the data on the permutation tests for allele frequencies for Table S5
gene_analysis/allele_frequency_permutations.RData : gene_analysis/GENE_OVERLAP_ANALYSIS \
	gene_analysis/allele_frequency_permutations.R
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

# Generating the gene ontology enrichment data for tables S6, S7 and S8
# Outputs the following files :
# bp_over_summary.RData
# bp_under_summary.RData
# pfam_over_summary.RData
gene_analysis/GO_ANALYSIS : gene_analysis/GENE_OVERLAP_ANALYSIS \
	gene_analysis/go_analysis.R \
	gene_analysis/soybase_genome_annotation_v4.0_04-20-2021.txt
	cd gene_analysis/ ; $(R_RUN_COMMAND) go_analysis.R; touch GO_ANALYSIS

# Generating the permutation test data for figure 5
gene_analysis/permutation_all_100kb.RData : gene_analysis/gprop_permutations.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	sv_genotyping/combined_svs/combined_paragraph_filtered.vcf \
	gene_analysis/GENE_OVERLAP_ANALYSIS
	cd gene_analysis ; $(R_RUN_COMMAND) gprop_permutations.R


##### TRANSPOSABLE ELEMENT ANALYSES

# Extracting a VCF of the sequences to query against the soybean TE database with BLASTN
te_analysis/query_all.vcf : sv_genotyping/combined_svs/combined_paragraph_filtered.vcf \
	te_analysis/extract_query_vcf.sh
	cd te_analysis ; ./extract_query_vcf.sh

# Querying the SV sequences against the SoyTEdb database retrieved from SoyBase
te_analysis/blast_svs.txt : te_analysis/query_all.vcf \
	te_analysis/blast_tes.sh \
	te_analysis/te_database/SoyBase_TE_Fasta.txt
	cd te_analysis ; ./blast_tes.sh $(BLASTN) $(MAKEBLASTDB)

# Identifying polymoprhic TEs within the combined Illumina/Oxford Nanopore SV dataset
te_analysis/polymorphic_tes.tsv : te_analysis/query_all.vcf \
	te_analysis/blast_svs.txt \
	te_analysis/te_blast_analysis.R \
	te_analysis/extract_te_metadata.sh \
	te_analysis/te_database/SoyBase_TE_Fasta.txt
	cd te_analysis ; $(R_RUN_COMMAND) te_blast_analysis.R

# Extracting VCF files corresponding to DNA TE SVs that have at least 3 matches in the dataset
te_analysis/DNA_TE_VCFS : te_analysis/query_all.vcf \
	te_analysis/polymorphic_tes.tsv \
	te_analysis/extract_dna_te_vcfs.R
	cd te_analysis ; $(R_RUN_COMMAND) extract_dna_te_vcfs.R ; touch DNA_TE_VCFS

# Assembling the sequences around DNA TE SVs and performing multiple alignments with MAFFT (ginsi command)
te_analysis/multiple_alignments/MULTIPLE_ALIGNMENTS : te_analysis/DNA_TE_VCFS \
	nanopore_data/NANOPORE_ALIGNMENT \
	te_analysis/multiple_alignments/make_symbolic_links.sh \
	te_analysis/multiple_alignments/assemble_align_all.R \
	te_analysis/multiple_alignments/age_realign.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	cd te_analysis/multiple_alignments ; ./make_symbolic_links.sh ; \
		$(R_RUN_COMMAND) assemble_align_all.R $(SAMTOOLS) $(MINIMAP2) $(WTDBG2) $(WTPOA_CNS) $(GINSI) ; touch MULTIPLE_ALIGNMENTS

# Filtering the multiple alignments, keeping only those with a conspicuous presence/absence polymorphism
# Also removing aligned sequences that have poor (< 0.75) identity to the reference sequence in their 500-nucleotide ends
te_analysis/multiple_alignments/FILTERED_ALIGNMENTS : te_analysis/multiple_alignments/MULTIPLE_ALIGNMENTS \
	te_analysis/multiple_alignments/filter_alignments.R \
	te_analysis/multiple_alignments/mafft_metadata.txt
	cd te_analysis/multiple_alignments ; $(R_RUN_COMMAND) filter_alignments.R ; touch FILTERED_ALIGNMENTS

# Detecting TIR and TSD sequences from the multiple alignments and analyzing the similarity of terminal repeats
te_analysis/multiple_alignments/TIR_TSD_ANALYSIS : te_analysis/multiple_alignments/FILTERED_ALIGNMENTS \
	te_analysis/multiple_alignments/mafft_metadata.txt \
	te_analysis/multiple_alignments/extract_te_sequences.R \
	te_analysis/multiple_alignments/grf_analysis.R \
	te_analysis/multiple_alignments/analyse_tir_sequences.R \
	te_analysis/multiple_alignments/filtered_alignments/correct_tsd_tir_sequences.txt
	cd te_analysis/multiple_alignments ; $(R_RUN_COMMAND) extract_te_sequences.R ; \
		$(R_RUN_COMMAND) grf_analysis.R $(GRF) ; \
		$(R_RUN_COMMAND) analyse_tir_sequences.R ; \
		touch TIR_TSD_ANALYSIS

te_analysis/multiple_alignments/Gm04_2257090_INS_480_analysis/STOWAWAY_MITE_ANALYSIS : structure_analysis/platypus_snps.vcf \
	te_analysis/query_all.vcf \
	te_analysis/multiple_alignments/Gm04_2257090_INS_480_analysis/stowaway_mite_analysis.R
	cd te_analysis/multiple_alignments/Gm04_2257090_INS_480_analysis ; \
		$(R_RUN_COMMAND) stowaway_mite_analysis.R $(BCFTOOLS) $(PLINK) $(BGZIP) $(TABIX) ; \
		touch STOWAWAY_MITE_ANALYSIS

