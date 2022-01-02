#!/prg/R/4.0/bin/Rscript

# Figure s7 shows the results of the subsampling analysis on insertion sensitivity and precision

# Loading the required packages
library(ggplot2)
library(grid)

# Loading the benchmark plotting data.frame that contains the processed data for plotting
# DEPENDENCY : nanopore_sv_calling/subsampling_analysis/benchmark_data.RData
load("../nanopore_sv_calling/subsampling_analysis/benchmark_data.RData")

# A function that takes the plotting dataframe, an SV type, and a measure (sensitivity or precision)
# and plots the results
plot_downsampling <- function(x, svtype, measure) {

	stopifnot(svtype %in% c("DEL", "INS", "INV", "DUP"))
	stopifnot(measure %in% c("sensitivity", "precision"))

	# Removing underscores in the cultivar names
	x$cultivar <- gsub("_", " ", x$cultivar)

	x <- x[x$size_class != "all" & x$svtype == svtype, ]

	if(measure == "sensitivity") {
		median_col = "s_median"
		min_col = "s_min"
		max_col = "s_max"
	} else {
		median_col = "p_median"
		min_col = "p_min"
		max_col = "p_max"
	
	}

	svtype_mapping <- c("DEL" = "deletions",
			    "INS" = "insertions",
			    "INV" = "inversions",
			    "DUP" = "duplications")

	size_class_levels <- c(paste0("[50-100 bp[ ", svtype_mapping[svtype]),
			       paste0("[100-1,000 bp[ ", svtype_mapping[svtype]),
			       paste0("[1,000-10,000 bp[ ", svtype_mapping[svtype]),
			       paste0("[10,000+ bp[ ", svtype_mapping[svtype]))

	size_class_mapping <- c("[50-100[" = size_class_levels[1],
				"[100-1000[" = size_class_levels[2],
				"[1000-10000[" = size_class_levels[3],
				"[10000+[" = size_class_levels[4])

	x$size_class <- factor(size_class_mapping[x$size_class], levels = size_class_levels)

	baseplot <- 
		ggplot(x[x$frac != 1, ], aes_string(x = "depth", y = median_col, ymin = min_col, ymax = max_col, color = "cultivar")) +
		geom_pointrange() +
		geom_point(data = x[x$frac == 1, ], shape = 8) +
		facet_wrap(~size_class, ncol = 1) +
		scale_y_continuous(name = ifelse(measure == "sensitivity", "Sensitivity", "Precision"),
				   limits = c(0, 1)) +
		xlab("Sequencing depth (X)") +
		theme_bw() +
		theme(panel.grid.minor = element_blank())

	main_plot <- baseplot + guides(color =  "none")
	legend_plot <- baseplot + theme(legend.position = "top", legend.direction = "horizontal")

	return(list(main_plot, legend_plot))
}

# Another function that takes a sensitivity plot, a precision plot, and a plot with a legend,
#  and wraps them up in a single plot
combine_plots <- function(sensitivity, precision, legend_plot) {
	grid.newpage()

	# Printing the plot with the legend first and then hiding everything but the legend
	print(legend_plot)
	pushViewport(viewport(y = 0, height = 0.95, just = "bottom", layout = grid.layout(nrow = 1, ncol = 2)))
	grid.rect(gp = gpar(col = "transparent", fill = "white"))

	# Printing the sensitivity plot
	pushViewport(viewport(layout.pos.col = 1))
	grid.text("A", x = 0.08, y = 0.96, gp = gpar(fontsize = 24))
	print(sensitivity, vp = viewport(x = 0.12, width = 0.88, just = "left"))
	popViewport()

	# Printing the precision plot
	pushViewport(viewport(layout.pos.col = 2))
	grid.text("B", x = 0.08, y = 0.96, gp = gpar(fontsize = 24))
	print(precision, vp = viewport(x = 0.12, width = 0.88, just = "left"))
	popViewport()
}

# A plot of deletions combining all the information
sensitivity_plot <- plot_downsampling(benchmark_plotting, "INS", "sensitivity")
precision_plot <- plot_downsampling(benchmark_plotting, "INS", "precision")

# theme parameters common to both subplots
common_theme <- theme(text = element_text(size = 14),
		      panel.grid.minor = element_blank())

# Saving the figure as Figure s7
png("figure_s7.png", width = 7, height = 8, units = "in", res = 500)
combine_plots(sensitivity_plot[[1]] + common_theme,
	      precision_plot[[1]] + common_theme, 
	      sensitivity_plot[[2]] + common_theme)
dev.off()

