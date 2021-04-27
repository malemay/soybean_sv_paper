#!/prg/R/4.0/bin/Rscript

# Code to create Figure S11 of the manuscript
# This figure is a scatter plot of the precision of insertion genotyping
#  as a function of Oxford Nanopore sequencing N50
# The objective is to show that samples that were sequenced with longer reads
#  had more power for insertion discovery, thus explaining the poor precision
#  of insertion genotyping for some samples
# Panel A will show results for insertions in the range 1,000-10,000, while
#  Panel B will show results for insertions > 10,000

# Loading the ggplot2 and grid packages
library(ggplot2)
library(grid)

# Loading the N50 data
# DEPENDENCY : N50_stats.txt
N50_stats <- read.table("/home/malem420/sv_manuscript/nanoplot_stats/N50_stats.txt", header = TRUE, stringsAsFactors = FALSE)

# Loading the genotyping precision and sensitivity rates
# DEPENDENCY : sveval_nogeno_rates.RData
load("/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_nanopore/genotyping/svmerged/paragraph/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData")

# We also load the correspondence between line names and their CAD IDs
# DEPENDENCY : line_ids.txt
line_ids <- read.table("/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_genotyping/svmerged_variants/paragraph/sveval_benchmarks/line_ids.txt", 
		       header = TRUE, stringsAsFactors = FALSE)

# Now we can extract and link the data through the IDs
line_ids$N50 <- N50_stats[match(line_ids$name, N50_stats$sample), "N50"]

# Now extracting the precision rates for ranges 1,000-10,000 and 10,000+
insertion_data <- sveval_nogeno_rates$INS
insertion_1000_data <- insertion_data[insertion_data$size_class == "[1000-10000[" & insertion_data$threshold == 2, ]
insertion_1000_precision <- insertion_1000_data$precision
names(insertion_1000_precision) <- insertion_1000_data$cultivar

insertion_10000_data <- insertion_data[insertion_data$size_class == "[10000+[" & insertion_data$threshold == 2, ]
insertion_10000_precision <- insertion_10000_data$precision
names(insertion_10000_precision) <- insertion_10000_data$cultivar

line_ids$insertions_1000 <- insertion_1000_precision[line_ids$id]
line_ids$insertions_10000 <- insertion_10000_precision[line_ids$id]

# Creating panel A (insertions in the range 1,000-10,000 bp)
panelA <- ggplot(line_ids, aes(x = N50 / 1000, y = insertions_1000)) +
	geom_point(size = 1) +
	scale_x_continuous(name = "Oxford Nanopore sequencing N50 (kb)") +
	scale_y_continuous(name = "Precision of insertion genotyping") +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
	      text = element_text(size = 10))

# And now creating panel B (insertions larger than 10,000 bp)
panelB <- ggplot(line_ids, aes(x = N50 / 1000, y = insertions_10000)) +
	geom_point(size = 1) +
	scale_x_continuous(name = "Oxford Nanopore sequencing N50 (kb)") +
	scale_y_continuous(name = "Precision of insertion genotyping") +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
	      text = element_text(size = 10))

# Saving to the png format
png("figure_s11.png", width = 6, height = 3, units = "in", res = 500)
pushViewport(viewport(layout = grid.layout(1, 2)))

# Printing panel A in the left column viewport
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(panelA, vp = viewport(x = 0.03, width = 0.97, just = "left"))
grid.text("A", x = 0.04, y = 0.96, gp = gpar(fontsize = 18))
popViewport()

# Printing panel B in the right column viewport
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(panelB, vp = viewport(x = 0.03, width = 0.97, just = "left"))
grid.text("B", x = 0.04, y = 0.96, gp = gpar(fontsize = 18))
popViewport()

dev.off()

