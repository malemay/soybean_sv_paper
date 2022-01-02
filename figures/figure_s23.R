#!/prg/R/4.0/bin/Rscript

# Figure s23 shows the phylogenetic trees obtained from SNPs and SVs

# Loading the ape package
library(ape)

# Loading the tree data
load("~/analyse_nanopore/manuscript_revision/population_analysis/tree/tree_data.RData")

# Saving the figure as Figure sA
png("figure_s23.png", width = 9.5, height = 6, units = "in", res = 500)

opar <- par(mar = c(0, 0, 0, 0), mfrow = c(1, 2))

plot(rotate(rotate(sv_tree, 124), 160), type = "fan", tip.color = pop_colors[sv_tree$tip.label], font = 1, cex = 0.45)
nodelabels(sv_support, cex = 0.4, bg = "white")
text(x = -0.28, y = 0.28, labels = "A", cex = 2.5)

plot(snp_tree, tip.color = pop_colors[snp_tree$tip.label], type = "fan", font = 1, cex = 0.45)
nodelabels(snp_support, cex = 0.4, bg = "white")
text(x = -0.28, y = 0.32, labels = "B", cex = 2.5)

par(opar)

dev.off()

