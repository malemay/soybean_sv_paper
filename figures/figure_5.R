#!/prg/R/4.0/bin/Rscript

# Figure 5 shows the allele frequencies of various SVs depending
# on the gene features (intergenic, 5 kb upstream, gene, exon)
# that they overlap

# Loading the ggplot2 library
library(ggplot2)

# First reading the dataset of overlapping features and frequencies
# DEPENDENCY : overlap_data.txt
overlap_data <- read.table("../gene_analysis/overlap_data.txt",
			   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

overlap_data$overlap <- factor(overlap_data$overlap, levels = c("cds", "gene", "upstream5kb", "intergenic"))

# Loading the data on the randomly shuffled SVs and gene overlaps
# DEPENDENCY : ../gene_analysis/permutation_all_100kb.RData
load("../gene_analysis/permutation_all_100kb.RData")

# Creating the plot
figure_5 <- ggplot(overlap_data, aes(x = af, color = overlap)) +
	geom_density() +
	facet_wrap(~svtype, ncol = 1) +
	scale_x_continuous(name = "Structural variant frequency") +
	scale_y_continuous(name = "Density") +
	scale_color_discrete(name = "Overlapped feature") +
	theme_bw() +
	theme(text = element_text(size = 14),
	      legend.position = "top",
	      legend.direction = "horizontal",
	      legend.key.size = unit(0.03, "npc"))

# Saving as a png file
png("figure_5.png", width = 6, height = 6, units = "in", res = 500)
print(figure_5)
dev.off()

