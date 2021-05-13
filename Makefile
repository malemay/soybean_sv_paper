
# Creating some variables for executables
# Using R 3.5.0 for figures because of a problem with resolution when using R 4.0
R_FIG_COMMAND = /usr/bin/R CMD BATCH --no-save --no-restore --no-site-file --no-environ
R_RUN_COMMAND = /prg/R/4.0/bin/R CMD BATCH --no-save --no-restore --no-site-file --no-environ
CIRCOS = /home/malem420/programs/circos-0.69-9/bin/circos

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

# --- The next section prepares the Illumina SV benchmarks from the Paragraph vcfs
#
ILLUMINA_BENCHMARK_VCFS = $(shell tail -n+2 utilities/line_ids.txt | cut -f2 | xargs -I {} echo sv_genotyping/illumina_svs/{}_results/genotypes.vcf.gz)
NANOPORE_NORMALIZED_SVS = $(shell tail -n+2 utilities/line_ids.txt | cut -f1 | xargs -I {} echo nanopore_sv_calling/{}_normalized_ids.vcf)

# Benchmark of Illumina SVs in non-repeat regions
sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData: \
	sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_benchmark.R \
	utilities/line_ids.txt \
	$(ILLUMINA_BENCHMARK_VCFS) \
	$(NANOPORE_NORMALIZED_SVS) \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/repeat_regions/non_repeated_regions.bed \
	scripts/extract_rates.R \
	scripts/read_filter_vcf.R
	cd sv_genotyping/illumina_svs/sveval_benchmarks ; $(R_RUN_COMMAND) norepeat_benchmark.R


