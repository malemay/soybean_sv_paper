#!/prg/R/4.0/bin/Rscript

# Figure S13 shows the benchmarking of duplications and inversions discovered using Oxford Nanopore sequencing
#  and genotyped using the Illumina sequencing data
# We only need a single panel to show inversions and another to show duplications because there are not many

# Loading the ggplot2 and grid packages
library(ggplot2)
library(grid)

# Loading the data used for plotting
# DEPENDENCY : sveval_nogeno_rates.RData
load("/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_nanopore/genotyping/svmerged/paragraph/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData")

# Also loading a script that will be used to prepare the data for plotting
# DEPENDENCY : make_plot_data.R
source("/home/malem420/scripts/make_plot_data.R")

# Preparing the data for plotting
dup_plot_data <- make_plot_data(sveval_nogeno_rates, "DUP")
inv_plot_data <- make_plot_data(sveval_nogeno_rates, "INV")

# We keep only the summaries over all size classes
dup_plot_data <- dup_plot_data[dup_plot_data$size_class == "all", ]
inv_plot_data <- inv_plot_data[inv_plot_data$size_class == "all", ]

# Defining a common theme for both plots
common_theme <- theme_bw() + 
	theme(plot.title = element_text(size = 12),
	      panel.grid.minor = element_blank(),
	      axis.title = element_blank())

# Defining common x- and y-axes
x_axis <- scale_x_continuous(name = "Sensitivity",
			     limits = c(0, 0.42), 
			     breaks = seq(0, 0.4, 0.1))

y_axis <- scale_y_continuous(name = "Precision",
			     limits = c(0, 1),
			     breaks = seq(0, 1, 0.2))

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
	common_theme +
	theme(axis.text.x = element_blank(),
	      axis.ticks.x = element_blank())

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


# Saving to disk as "figure_s13.png"
png("figure_s13.png", width = 3, height = 6, units = "in", res = 500)
grid.newpage()
# Locating the subplots in the figure
pushViewport(viewport(x = 0.05, y = 0.03, just = c("left", "bottom"), width = 0.95, height = 0.97))
pushViewport(viewport(layout = grid.layout(2, 1)))
print(duplications_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(inversions_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

# Adding axis text
popViewport(2)
pushViewport(viewport(x = 0, y = 0.03, just = c("left", "bottom"), width = 0.05, height = 0.97))
grid.text("Precision", rot = 90)
popViewport()
pushViewport(viewport(x = 0.05, y = 0, just = c("left", "bottom"), width = 0.95, height = 0.03))
grid.text("Sensitivity")

dev.off()

