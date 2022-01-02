#!/usr/bin/Rscript

# Figure 1 shows the benchmarking of SVs discovered using Illumina sequencing
# For deletions, 4 panels are needed (50-100, 100-1000, 1000-10000 and 10000+)
# For insertions, only 2 panels are needed because almost no SVs > 1000 bp were found

# Loading the ggplot2 and grid packages
library(ggplot2)
library(grid)

# Loading the data used for plotting
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData
load("../sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData")

# Also loading a script that will be used to prepare the data for plotting
# DEPENDENCY : scripts/make_plot_data.R
source("../scripts/make_plot_data.R")

# Preparing the data for plotting
del_plot_data <- make_plot_data(sveval_nogeno_rates, "DEL")
ins_plot_data <- make_plot_data(sveval_nogeno_rates, "INS")

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
	theme(text = element_text(size = 10),
	      strip.background = element_blank(),
	      panel.grid.minor = element_blank(),
	      legend.text = element_text(size = 6),
	      legend.key.height = unit(0.02, "npc"),
	      legend.position = "top",
	      legend.direction = "horizontal")

# Defining a common y-axis
y_axis <- scale_y_continuous(name = "Precision",
			     limits = c(0.4, 1),
			     expand = c(-0.02, 0.02),
			     breaks = seq(0.4, 1, 0.2))

# Preparing the plot for deletions
deletions_plot <- 
	ggplot(del_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.2) +
	geom_point(mapping = aes(color = cultivar), size = 0.3) +
	geom_point(data = del_plot_data[del_plot_data$threshold == 2, ], aes(color = cultivar), shape = 8, size = 1.5) +
	# Adding labels with the number of supporting reads for a single sample in one panel
	geom_text(data = del_plot_data[del_plot_data$size_class == "[50-100[" & del_plot_data$cultivar == "CAD1052" &
		  del_plot_data$threshold %in% c(2, 6, 10, 14, 16, 18, 20, 22, 24, 26, 29), ],
		  aes(label = threshold, color = cultivar, y = precision + 0.04), size = 2) +
	facet_wrap(~size_class,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ deletions",
					 "[100-1000[" = "[100-1,000 bp[ deletions",
					 "[1000-10000[" = "[1,000-10,000 bp[ deletions",
					 "[10000+[" = "[10,000+ bp[ deletions"))) +
	scale_x_continuous(name = "Sensitivity",
			   limits = c(-0.02, 0.82), 
			   expand = c(-0.02, 0.02),
			   breaks = seq(0, 0.8, 0.2)) +
	y_axis +
	guides(color = FALSE) +
	common_theme +
	theme(panel.spacing.y = unit(0.03, "npc"))

# Now preparing the plot for insertions
insertions_plot <- 
	ggplot(ins_plot_data, aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = cultivar), size = 0.2) +
	geom_point(mapping = aes(color = cultivar), size = 0.5) +
	geom_point(data = ins_plot_data[ins_plot_data$threshold == 2, ], aes(color = cultivar), shape = 8, size = 2) +
	facet_wrap(~size_class,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ insertions",
					 "[100-1000[" = "[100-1,000 bp[ insertions"))) +
	scale_x_continuous(name = "Sensitivity") +
	y_axis +
	guides(color = FALSE) +
	common_theme


# Saving to disk as "Figure_1.png"
# OUTPUT : figures/figure_1.png
png("figure_1.png", width = 10, height = 16, units = "cm", res = 800)
grid.newpage()
# Locating the subplots in the figure
pushViewport(viewport(x = 0.03, width = 0.97, just = "left"))
pushViewport(viewport(layout = grid.layout(2, 1, heights = c(unit(65, "null"), unit(35, "null")))))
del_vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
print(deletions_plot, vp = del_vp)
ins_vp <- viewport(layout.pos.row = 2, layout.pos.col = 1)
print(insertions_plot, vp = ins_vp)
# Adding the A and B panel labels
grid.text("A", x = 0, y = 0.97, gp = gpar(fontsize = 16), vp = del_vp)
grid.text("B", x = 0, y = 0.97, gp = gpar(fontsize = 16), vp = ins_vp)
dev.off()

