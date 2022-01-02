#!/prg/R/4.0/bin/Rscript

# Figure S12 shows the benchmarking of SVs discovered using Illumina sequencing
# For deletions, 4 panels are needed (50-100, 100-1000, 1000-10000 and 10000+)
# For insertions, only 2 panels are needed because almost no SVs > 1000 bp were found

# This figure is similar to figures 1 and S4 except that here, the quality measure that
# is used for the precision-recall curve is the number of independent programs that called
# the variant, while the DP (number of supporting reads) threshold is maintained to 2 and
# the AC_Hom (the number of ALT alleles within homozygous ALT calls) threshold is maintained
# to 4.

# Loading the ggplot2 and grid packages
library(ggplot2)
library(grid)

# Loading the data used for plotting
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/sveval_ncallers_rates.RData
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/NCALLERS_ILLUMINA_BENCHMARK
load("../sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/sveval_ncallers_rates.RData")

# Also loading a script that will be used to prepare the data for plotting
# DEPENDENCY : scripts/make_plot_data.R
source("../scripts/make_plot_data.R")

# Preparing the data for plotting
del_plot_data <- make_plot_data(sveval_ncallers_rates, "DEL")
ins_plot_data <- make_plot_data(sveval_ncallers_rates, "INS")

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
x_axis <- scale_x_continuous(name = "Sensitivity")
	#		     limits = c(0, 0.8), 
	#		     breaks = seq(0, 0.8, 0.2))

y_axis <- scale_y_continuous(name = "Precision")
			     #limits = c(0.6, 1),
			     #expand = c(-0.02, 0.02),
			     #breaks = seq(0.6, 1, 0.1))

# Preparing the plot for deletions
deletions_plot <- 
	ggplot(del_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.2) +
	geom_point(mapping = aes(color = cultivar), size = 0.5) +
	geom_point(data = del_plot_data[del_plot_data$threshold == 1, ], aes(color = cultivar), shape = 8, size = 2) +
	facet_wrap(~size_class,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ deletions",
					 "[100-1000[" = "[100-1,000 bp[ deletions",
					 "[1000-10000[" = "[1,000-10,000 bp[ deletions",
					 "[10000+[" = "[10,000+ bp[ deletions"))) +

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
	geom_point(data = ins_plot_data[ins_plot_data$threshold == 1, ], aes(color = cultivar), shape = 8, size = 2) +
	facet_wrap(~size_class,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ insertions",
					 "[100-1000[" = "[100-1,000 bp[ insertions"))) +
	scale_x_continuous(name = "Sensitivity",
			   limits = c(0, 0.42),
			   breaks = seq(0, 0.4, 0.1)) +
	y_axis +
	guides(color = FALSE) +
	common_theme

# Saving as a png file
# OUTPUT : figures/figure_s12.png
png("figure_s12.png", width = 6, height = 9, units = "in", res = 500)
grid.newpage()
# Locating the subplots in the figure
pushViewport(viewport(x = 0.05, just = "left", width = 0.95))
pushViewport(viewport(layout = grid.layout(2, 1, heights = c(unit(65, "null"), unit(35, "null")))))
del_vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
print(deletions_plot, vp = del_vp)
ins_vp <- viewport(layout.pos.row = 2, layout.pos.col = 1)
print(insertions_plot, vp = ins_vp)
# Adding the A and B panel labels
grid.text("A", x = 0, y = 0.97, gp = gpar(fontsize = 24), vp = del_vp)
grid.text("B", x = 0, y = 0.97, gp = gpar(fontsize = 24), vp = ins_vp)
dev.off()

