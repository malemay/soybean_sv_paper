#!/usr/bin/Rscript

# Figure 2 shows the benchmarking of SVs discovered using Oxford Nanopore sequencing
#  but genotyped using Illumina sequencing data
# Four (4) panels are needed for both deletions and insertions because we want to show
#  the results for all size classes for both SV types

# Loading the ggplot2 and grid packages
library(ggplot2)
library(grid)

# DEPENDENCY : sv_genotyping/nanopore_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData
# Loading the data used for plotting
load("../sv_genotyping/nanopore_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData")

# Also loading a script that will be used to prepare the data for plotting
# DEPENDENCY : scripts/make_plot_data.R
source("../scripts/make_plot_data.R")

# Preparing the data for plotting
del_plot_data <- make_plot_data(sveval_nogeno_rates, "DEL")
ins_plot_data <- make_plot_data(sveval_nogeno_rates, "INS")

# We remove the "[30-50[" and "all" size classes
del_plot_data <- del_plot_data[!del_plot_data$size_class %in% c("[30-50[", "all"), ]
del_plot_data$size_class <- droplevels(del_plot_data$size_class)

ins_plot_data <- ins_plot_data[!ins_plot_data$size_class %in% c("[30-50[", "all"), ]
ins_plot_data$size_class <- droplevels(ins_plot_data$size_class)

# Defining a common theme for both plots
common_theme <- theme_bw() + 
	theme(strip.background = element_blank(),
	      panel.grid.minor = element_blank())

# Defining common x- and y-axes
x_axis <- scale_x_continuous(name = "Sensitivity",
			     limits = c(0, 1), 
			     breaks = seq(0, 1, 0.2))

y_axis <- scale_y_continuous(name = "Precision",
			     limits = c(0, 1.03),
			     breaks = seq(0, 1, 0.2))

# Preparing the plot for deletions
deletions_plot <- 
	ggplot(del_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.2) +
	geom_point(mapping = aes(color = cultivar), size = 0.5) +
	geom_point(data = del_plot_data[del_plot_data$threshold == 2, ], aes(color = cultivar), shape = 8, size = 2) +
	# Adding labels with the number of supporting reads for a single sample in one panel
	geom_text(data = del_plot_data[del_plot_data$size_class == "[50-100[" & del_plot_data$cultivar == "CAD1052" &
		  del_plot_data$threshold %in% c(2, 6, 8, 12, 16, 19, 23, 30), ],
		  aes(label = threshold, color = cultivar, y = precision + 0.06), size = 2) +
	facet_wrap(~size_class, ncol = 2,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ deletions",
					 "[100-1000[" = "[100-1,000 bp[ deletions",
					 "[1000-10000[" = "[1,000-10,000 bp[ deletions",
					 "[10000+[" = "[10,000+ bp[ deletions"))) +
	x_axis +
	y_axis +
	guides(color = FALSE) +
	common_theme

# Now preparing the plot for insertions
insertions_plot <- 
	ggplot(ins_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.2) +
	geom_point(mapping = aes(color = cultivar), size = 0.5) +
	geom_point(data = ins_plot_data[ins_plot_data$threshold == 2, ], aes(color = cultivar), shape = 8, size = 2) +
	facet_wrap(~size_class,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ insertions",
					 "[100-1000[" = "[100-1,000 bp[ insertions",
					 "[1000-10000[" = "[1,000-10,000 bp[ insertions",
					 "[10000+[" = "[10,000+ bp[ insertions"))) +
	x_axis +
	y_axis +
	guides(color = FALSE) +
	common_theme


# Saving to disk as "figure_2.png"
# OUTPUT : figures/figure_2.png
png("figure_2.png", width = 10, height = 20, units = "cm", res = 800)
grid.newpage()
# Locating the subplots in the figure, leaving some space for the "A" and "B" plot labels
pushViewport(viewport(layout = grid.layout(2, 1)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(deletions_plot, vp = viewport(x = 0.03, width = 0.97, just = "left"))
grid.text("A", x = 0.04, y = 0.95, gp = gpar(fontsize = 16))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
print(insertions_plot, vp = viewport(x = 0.03, width = 0.97, just = "left"))
grid.text("B", x = 0.04, y = 0.95, gp = gpar(fontsize = 16))
popViewport()

dev.off()

