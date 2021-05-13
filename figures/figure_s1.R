#!/prg/R/4.0/bin/Rscript

# Figure S1 shows the benchmarking of SVs discovered using Illumina sequencing
# For deletions, 4 panels are needed (50-100, 100-1000, 1000-10000 and 10000+)
# For insertions, only 2 panels are needed because almost no SVs > 1000 bp were found

# The difference with Figure 1 is that we look at non-repeated regions
# Analysis was restricted to those regions using the bed.regions argument
# of sveval::svevalOl

# Loading the ggplot2 and grid packages
library(ggplot2)
library(grid)

# Loading the data used for plotting
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData
load("../sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData/sveval_norepeat_rates.RData")

# Also loading a script that will be used to prepare the data for plotting
# DEPENDENCY : scripts/make_plot_data.R
source("../scripts/make_plot_data.R")

# Preparing the data for plotting
del_plot_data <- make_plot_data(sveval_norepeat_rates, "DEL")
ins_plot_data <- make_plot_data(sveval_norepeat_rates, "INS")

# Adding the F-score to the data.frame
del_plot_data$Fscore <- 2 * (del_plot_data$precision * del_plot_data$sensitivity) / (del_plot_data$precision + del_plot_data$sensitivity)
ins_plot_data$Fscore <- 2 * (ins_plot_data$precision * ins_plot_data$sensitivity) / (ins_plot_data$precision + ins_plot_data$sensitivity)

# For deletions, F-score appears to be maximized at a threshold of 2 in most cases
# For deletions, F-score appears to be maximized at a threshold of 2 in most cases;
#  however, I will use a filtering threshold of 2 for both the sake of consistency
#  across variant types and because a threshold of 1 read sounds extremely lenient

# For deletions, we remove the "[30-50[" and "all" columns
del_plot_data <- del_plot_data[!del_plot_data$size_class %in% c("[30-50[", "all"), ]
del_plot_data$size_class <- droplevels(del_plot_data$size_class)

# For insertions, we only keep the "[50-100[" and [100-1000[" columns
ins_plot_data <- ins_plot_data[ins_plot_data$size_class %in% c("[50-100[", "[100-1000["), ]
ins_plot_data$size_class <- droplevels(ins_plot_data$size_class)


# Defining a common theme for both plots
common_theme <- theme_bw() + 
	theme(panel.grid.minor = element_blank(),
	      text = element_text(size = 15))

# Defining common x- and y-axes
x_axis <- scale_x_continuous(name = "Sensitivity",
			     limits = c(-0.02, 0.85), 
			     breaks = seq(0, 0.8, 0.2))

y_axis <- scale_y_continuous(name = "Precision",
			     limits = c(0.4, 1),
			     breaks = seq(0.4, 1, 0.2))

# Preparing the plot for deletions
deletions_plot <- 
	ggplot(del_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.2) +
	geom_point(mapping = aes(color = cultivar), size = 0.5) +
	geom_point(data = del_plot_data[del_plot_data$threshold == 2, ], aes(color = cultivar), shape = 8, size = 2) +
	facet_wrap(~size_class,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ deletions",
					 "[100-1000[" = "[100-1,000 bp[ deletions",
					 "[1000-10000[" = "[1,000-10,000 bp[ deletions",
					 "[10000+[" = "[10,000+ bp deletions"))) +

	x_axis +
	y_axis +
	guides(color = FALSE) +
	common_theme +
	theme( panel.spacing.y = unit(0.03, "npc"))

# Now preparing the plot for insertions
insertions_plot <- 
	ggplot(ins_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.2) +
	geom_point(mapping = aes(color = cultivar), size = 0.5) +
	geom_point(data = ins_plot_data[ins_plot_data$threshold == 2, ], aes(color = cultivar), shape = 8, size = 2) +
	facet_wrap(~size_class,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ insertions",
					 "[100-1000[" = "[100-1,000[ bp insertions"))) +
	x_axis +
	y_axis +
	guides(color = FALSE) +
	common_theme

# Saving in png format
# OUTPUT : figures/figure_s1.png
png("figure_s1.png", width = 6, height = 9, units = "in", res = 500)
grid.newpage()
# Locating the subplots in the figure
pushViewport(viewport(x = 0.03, just = "left", width = 0.97))
pushViewport(viewport(layout = grid.layout(2, 1, heights = c(unit(65, "null"), unit(35, "null")))))
del_vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
print(deletions_plot, vp = del_vp)
ins_vp <- viewport(layout.pos.row = 2, layout.pos.col = 1)
print(insertions_plot, vp = ins_vp)
# Adding the panel labels
grid.text("A", x = 0, y = 0.97, gp = gpar(fontsize = 24), vp = del_vp)
grid.text("B", x = 0, y = 0.97, gp = gpar(fontsize = 24), vp = ins_vp)
dev.off()

