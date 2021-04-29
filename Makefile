
# Creating some variables for executables
# Using R 3.5.0 for figures because of a problem with resolution when using R 4.0
R_FIG_COMMAND = /usr/bin/R CMD BATCH --no-save --no-restore --no-site-file --no-environ
R_RUN_COMMAND = /prg/R/4.0/bin/R CMD BATCH --no-save --no-restore --no-site-file --no-environ

# Creating some variables for more readable coding
SUPFIGURES = figures/figure_s1.png figures/figure_s2.png figures/figure_s3.png figures/figure_s4.png figures/figure_s5.png figures/figure_s6.png figures/figure_s7.png figures/figure_s8.png figures/figure_s9.png figures/figure_s10.png figures/figure_s11.png figures/figure_s12.png figures/figure_s13.png figures/figure_s14.png figures/figure_s15.png figures/figure_s16.png figures/figure_s17.png figures/figure_s18.png figures/figure_s19.png
SUPTABLES = tables/table_s1.csv tables/table_s2.csv tables/table_s3.csv tables/table_s4.csv
FIGURES = figures/figure_1.png figures/figure_2.png figures/figure_3_circos/figure_3.png figures/figure_4.png figures/figure_5.png figures/figure_6.png
TABLES = tables/table_1.png tables/table_2.png tables/table_3.png
CIRCOS = /home/malem420/programs/circos-0.69-9/bin/circos

# --- This target prepares all the figures, tables, and supplemental data
all: Supplemental_Data.pdf $(TABLES) $(FIGURES)

# --- The following section prepares de Supplemental Data file
Supplemental_Data.pdf : Supplemental_Data.tex references.bib genome_research.bst $(SUPFIGURES) $(SUPTABLES) 
	pdflatex Supplemental_Data.tex; bibtex Supplemental_Data; pdflatex Supplemental_Data.tex; pdflatex Supplemental_Data.tex

figures/figure_%.png : figures/figure_%.R
	cd figures; $(R_FIG_COMMAND) $(<F)

tables/table_s1.csv tables/table_s2.csv tables/table_s3.csv: tables/formatting_sup_tables.R
	cd tables; $(R_RUN_COMMAND) formatting_sup_tables.R

tables/table_s4.csv:
	cd tables; $(R_RUN_COMMAND) table_s4.R

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

figures/figure_3_circos/figure_3.png: figures/figure_3_circos/circos.conf # WARNING: the prerequisites here are incomplete
	cd figures/figure_3_circos/; $(CIRCOS)

figures/figure_4.png: figures/figure_4.R
	cd figures; $(R_FIG_COMMAND) figure_4.R

figures/figure_5.png: figures/figure_5.R
	cd figures; $(R_FIG_COMMAND) figure_5.R

figures/figure_6.png: figures/figure_6.R
	cd figures; $(R_FIG_COMMAND) figure_6.R

