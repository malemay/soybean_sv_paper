#!/usr/bin/Rscript

# Figure S17 shows the benchmarking of SVs discovered using Oxford Nanopore sequencing
#  but genotyped using Illumina sequencing data
# Four (4) panels are needed for both deletions and insertions because we want to show
#  the results for all size classes for both SV types
# The only difference with Figure 2 is that here we show the benchmarks obtained in non-repeat regions

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
del_plot_data <- make_plot_data(sveval_norepeat_rates, "DEL")
ins_plot_data <- make_plot_data(sveval_norepeat_rates, "INS")

# We remove the "[30-50[" and "all" size classes
del_plot_data <- del_plot_data[!del_plot_data$size_class %in% c("[30-50[", "all"), ]
del_plot_data$size_class <- droplevels(del_plot_data$size_class)

ins_plot_data <- ins_plot_data[!ins_plot_data$size_class %in% c("[30-50[", "all"), ]
ins_plot_data$size_class <- droplevels(ins_plot_data$size_class)

# Defining a common theme for both plots
common_theme <- theme_bw() + 
	theme(panel.grid.minor = element_blank(),
	      text = element_text(size = 16))

# Defining common x- and y-axes
x_axis <- scale_x_continuous(name = "Sensitivity",
			     limits = c(0, 1), 
			     breaks = seq(0, 1, 0.2))

y_axis <- scale_y_continuous(name = "Precision",
			     limits = c(0, 1),
			     breaks = seq(0, 1, 0.2))

# Preparing the plot for deletions
deletions_plot <- 
	ggplot(del_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.3) +
	geom_point(mapping = aes(color = cultivar), size = 1.5) +
	geom_point(data = del_plot_data[del_plot_data$threshold == 2, ], aes(color = cultivar), shape = 8, size = 2) +
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
	geom_line(mapping = aes(group = cultivar), size = 0.3) +
	geom_point(mapping = aes(color = cultivar), size = 1.5) +
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

# Saving to disk as "figure_s17.png"
# OUTPUT : figures/figure_s17.png
png("figure_s17.png", width = 12, height = 6, units = "in", res = 500)
grid.newpage()
# Locating the subplots in the figure, leaving some space for the "A" and "B" plot labels
pushViewport(viewport(layout = grid.layout(1, 2)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(deletions_plot, vp = viewport(x = 0.03, width = 0.97, just = "left"))
grid.text("A", x = 0.04, y = 0.97, gp = gpar(fontsize = 24))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(insertions_plot, vp = viewport(x = 0.03, width = 0.97, just = "left"))
grid.text("B", x = 0.04, y = 0.97, gp = gpar(fontsize = 24))
popViewport()

dev.off()

