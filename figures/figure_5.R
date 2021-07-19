#!/usr/bin/Rscript

# Figure 5 shows the allele frequencies of various SVs depending
# on the gene features (intergenic, 5 kb upstream, gene, exon)
# that they overlap

# Loading the ggplot2 library
library(ggplot2)
library(grid)

# First reading the dataset of overlapping features and frequencies
# DEPENDENCY : gene_analysis/overlap_data.txt
overlap_data <- read.table("../gene_analysis/overlap_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ordering the levels of the features
overlap_data$overlap <- factor(overlap_data$overlap, levels = c("cds", "gene", "upstream5kb", "intergenic"))

# Adding an allele frequency class to plot the number of SVs per bin
overlap_data$af_class <- cut(overlap_data$af, seq(0, 1, 0.1))


# Creating panel B of the plot
figure_5b <- ggplot(overlap_data[overlap_data$ af != 1, ], aes(x = af_class, fill = overlap)) +
	geom_bar(width = 0.6, position = position_dodge()) +
	facet_wrap(~svtype, ncol = 1,
		   labeller = labeller(svtype = c("DEL" = "Deletions", "INS" = "Insertions"))) +
	scale_x_discrete(name = "Allele frequency bin") +
	scale_y_log10(name = "Number of SVs") +
	scale_fill_discrete(name = "Overlapped feature") +
	theme_bw() +
	theme(text = element_text(size = 6),
	      axis.text = element_text(size = 4, color = "black"),
	      axis.ticks = element_line(size = 0.2),
	      strip.text = element_text(size = 6),
	      strip.background = element_blank(),
	      panel.spacing = unit(0.01, "npc"),
	      legend.position = "top",
	      legend.direction = "horizontal",
	      legend.box.margin = margin(0.005, 0.005, 0.005, 0.005, "npc"),
	      legend.box.spacing = unit(0, "npc"),
	      legend.key.size = unit(0.03, "npc"),
	      panel.grid.major = element_line(size = 0.2),
	      panel.grid.minor = element_blank())

# Loading the data on the randomly shuffled SVs and gene overlaps
# DEPENDENCY : gene_analysis/permutation_all_100kb.RData
load("../gene_analysis/permutation_all_100kb.RData")

# Formatting this data for plotting with ggplot2
del_distributions <- data.frame(proportion = as.numeric(permutation_all_100kb$del),
				feature = rep(colnames(permutation_all_100kb$del), each = nrow(permutation_all_100kb$del)),
				svtype = "DEL",
				stringsAsFactors = FALSE)

ins_distributions <- data.frame(proportion = as.numeric(permutation_all_100kb$ins),
				feature = rep(colnames(permutation_all_100kb$ins), each = nrow(permutation_all_100kb$ins)),
				svtype = "INS",
				stringsAsFactors = FALSE)

random_proportions <- rbind(del_distributions, ins_distributions)
random_proportions$feature <- factor(random_proportions$feature, levels = c("cds", "gene", "upstream5kb", "intergenic"))

# Creating a data.frame of observed values to add to the plot
gene_table <- as.matrix(table(overlap_data$overlap, overlap_data$svtype))
gene_table <- data.frame(DEL = gene_table[, "DEL"], INS = gene_table[, "INS"])
gene_table[, "DELprop"] <- gene_table[, "DEL"] / sum(gene_table[, "DEL"])
gene_table[, "INSprop"] <- gene_table[, "INS"] / sum(gene_table[, "INS"])
observed_proportions <- data.frame(proportion = c(gene_table$DELprop, gene_table$INSprop),
				   feature = rep(rownames(gene_table), 2),
				   svtype = rep(c("DEL", "INS"), each = nrow(gene_table)),
				   stringsAsFactors = FALSE)
observed_proportions$feature <- factor(observed_proportions$feature, levels = c("cds", "gene", "upstream5kb", "intergenic"))

# Creating panel A of the plot
figure_5a <- ggplot(random_proportions, aes(x = proportion)) +
	geom_histogram(fill = "white", color = "blue", bins = 60, size = 0.1) +
	geom_vline(data = observed_proportions, mapping = aes(xintercept = proportion), linetype = 3, color = "black", size = 0.15) +
	geom_blank(data = data.frame(proportion = c(0.09, 0.32), svtype = c("INS", "INS"), feature = c("gene", "upstream5kb"))) +
	facet_grid(svtype ~ feature, scales = "free_x",
		   labeller = labeller(svtype = c("DEL" = "Deletions", "INS" = "Insertions"))) +
	scale_x_continuous(name = "Proportion of SVs overlapping feature") +
	scale_y_continuous(name = "Count") +
	theme_bw() +
	theme(panel.grid = element_blank(),
	      axis.ticks = element_line(size = 0.2),
	      axis.text = element_text(size = 4, color = "black"),
	      text = element_text(size = 6), 
	      strip.text = element_text(size = 6,
					margin = margin(0, 2, 2, 2, "pt")),
	      strip.background = element_blank(),
	      panel.spacing = unit(0.01, "npc"))

# Arranging the two panels horizontally and saving as a png file
# OUTPUT : figures/figure_5.png
png("figure_5.png", width = 6, height = 3, units = "in", res = 500)
grid.newpage()

# Creating the two-panel layout
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
# Pushing the left viewport, leaving some space to the left for the panel label
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(figure_5a, vp = viewport(x = 0.02, width = 0.98, just = "left"))
grid.text("A", x = 0.05, y = 0.95)
popViewport()
# Doing the same thing with the right viewport
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(figure_5b, vp = viewport(x = 0.02, width = 0.98, just = "left"))
grid.text("B", x = 0.05, y = 0.95)
popViewport()

dev.off()

