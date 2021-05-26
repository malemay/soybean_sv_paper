#!/usr/bin/Rscript

# Figure S14 shows the benchmarking of duplications and inversions discovered using Oxford Nanopore sequencing
#  and genotyped using the Illumina sequencing data in non-repeat regions
# We only need a single panel to show inversions and another to show duplications because there are not many

# Loading the ggplot2 and grid packages
library(ggplot2)
library(grid)

# Loading the data used for plotting
# DEPENDENCY : sv_genotyping/nanopore_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData
load("../sv_genotyping/nanopore_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData")

# Also loading a script that will be used to prepare the data for plotting
# DEPENDENCY : scripts/make_plot_data.R
source("../scripts/make_plot_data.R")

# Preparing the data for plotting
dup_plot_data <- make_plot_data(sveval_norepeat_rates, "DUP")
inv_plot_data <- make_plot_data(sveval_norepeat_rates, "INV")

# We keep only the summaries over all size classes
dup_plot_data <- dup_plot_data[dup_plot_data$size_class == "all", ]
inv_plot_data <- inv_plot_data[inv_plot_data$size_class == "all", ]

# Defining a common theme for both plots
common_theme <- theme_bw() + 
	theme(plot.title = element_text(size = 12),
	      panel.grid.minor = element_blank())

# Defining common x- and y-axes
x_axis <- scale_x_continuous(name = "Sensitivity")

y_axis <- scale_y_continuous(name = "Precision")

# Preparing the plot for duplications
duplications_plot <- 
	ggplot(dup_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.2) +
	geom_point(mapping = aes(color = cultivar), size = 0.5) +
	geom_point(data = dup_plot_data[dup_plot_data$threshold == 2, ], aes(color = cultivar), shape = 8, size = 2) +
	x_axis +
	y_axis +
	guides(color = FALSE) +
	ggtitle("Duplications (all sizes)") +
	common_theme

# Now preparing the plot for inversions
inversions_plot <- 
	ggplot(inv_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.2) +
	geom_point(mapping = aes(color = cultivar), size = 0.5) +
	geom_point(data = inv_plot_data[inv_plot_data$threshold == 2, ], aes(color = cultivar), shape = 8, size = 2) +
	x_axis +
	y_axis +
	guides(color = FALSE) +
	ggtitle("Inversions (all sizes)") +
	common_theme


# Saving to disk as "figure_s14.png"
# OUTPUT : figures/figure_s14.png
png("figure_s14.png", width = 3, height = 6, units = "in", res = 500)
grid.newpage()
# Locating the subplots in the figure
pushViewport(viewport(x = 0.05, just = "left", width = 0.95))
pushViewport(viewport(layout = grid.layout(2, 1)))
dup_vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
print(duplications_plot, vp = dup_vp)
inv_vp <- viewport(layout.pos.row = 2, layout.pos.col = 1)
print(inversions_plot, vp = inv_vp)
# Add the panel labels A and B
grid.text("A", x = 0, y = 0.95, gp = gpar(fontsize = 20), vp = dup_vp)
grid.text("B", x = 0, y = 0.95, gp = gpar(fontsize = 20), vp = inv_vp)
dev.off()

