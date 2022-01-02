#!/usr/bin/Rscript

# Loading the required libraries
library(ggplot2)
library(grid)

# Loading the datasets
# DEPENDENCY : breakpoint_refinement_analysis/raw_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData
load("../breakpoint_refinement_analysis/raw_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData")
raw_variants_rates <- sveval_nogeno_rates

# DEPENDENCY : breakpoint_refinement_analysis/refined_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData
load("../breakpoint_refinement_analysis/refined_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData")
realigned_variants_rates <- sveval_nogeno_rates

rm(sveval_nogeno_rates)

# Sourcing the function used to prepare the data for plotting
# DEPENDENCY : scripts/make_plot_data.R
source("../scripts/make_plot_data.R")

# Plot-ready data.frames for deletions and insertions are prepared
raw_deletions <- make_plot_data(raw_variants_rates, "DEL")
raw_deletions$type <- "raw"
realigned_deletions <- make_plot_data(realigned_variants_rates, "DEL")
realigned_deletions$type <- "realigned"
deletions <- rbind(raw_deletions, realigned_deletions)
deletions$linegroup2 <- paste0(deletions$linegroup, deletions$type)

raw_insertions <- make_plot_data(raw_variants_rates, "INS")
raw_insertions$type <- "raw"
realigned_insertions <- make_plot_data(realigned_variants_rates, "INS")
realigned_insertions$type <- "realigned"
insertions <- rbind(raw_insertions, realigned_insertions)
insertions$linegroup2 <- paste0(insertions$linegroup, insertions$type)

# Removing OAC PETREL from the data and renaming the cultivars
deletions <- deletions[deletions$cultivar != "CAD1022", ]
insertions <- insertions[insertions$cultivar != "CAD1022", ]

line_names <- c("CAD1010" = "Maple Presto",
		"CAD1052" = "OAC Embro",
		"CAD1064" = "OAC Carman",
		"CAD1070" = "QS5091.50j")

deletions$cultivar <- line_names[deletions$cultivar]
insertions$cultivar <- line_names[insertions$cultivar]

# Removing the "[30-50[" and "all" size classes
deletions <- deletions[!deletions$size_class %in% c("[30-50[", "all"), ]
deletions$size_class <- droplevels(deletions$size_class)

insertions <- insertions[!insertions$size_class %in% c("[30-50[", "all"), ]
insertions$size_class <- droplevels(insertions$size_class)

# Defining a common theme for both plots
common_theme <- theme_bw() + 
	theme(panel.grid.minor = element_blank(),
	      text = element_text(size = 16),
	      legend.title = element_text(size = 12),
	      legend.text = element_text(size = 11),
	      legend.key.height = unit(0.02, "npc"),
	      legend.position = "top",
	      legend.direction = "horizontal")

# Defining common x- and y-axes
x_axis <- scale_x_continuous(name = "Sensitivity",
                             limits = c(0, 1), 
                             breaks = seq(0, 1, 0.2))

y_axis <- scale_y_continuous(name = "Precision",
                             limits = c(0, 1),
                             breaks = seq(0, 1, 0.2))

# Panel A shows the deletions
panelA <- 
	ggplot(deletions[deletions$pipeline == "paragraph", ], aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = linegroup2), size = 0.3) +
	geom_point(mapping = aes(color = type, shape = cultivar), size = 1.5) +
	facet_wrap(~size_class, ncol = 2,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ deletions",
					 "[100-1000[" = "[100-1,000 bp[ deletions",
					 "[1000-10000[" = "[1,000-10,000 bp[ deletions",
					 "[10000+[" = "[10,000+ bp[ deletions"))) +
	x_axis +
	y_axis +
	scale_color_discrete(name = "",
			     labels = c("realigned" = "Refined SVs", "raw" = "Raw SVs")) +
	guides(shape = FALSE,
	       color = guide_legend(override.aes = list(size = 3))) +
	common_theme

# Panel B shows the insertions
panelB <- 
	ggplot(insertions[insertions$pipeline == "paragraph", ], aes(x = sensitivity, y = precision)) +
	geom_line(mapping = aes(group = linegroup2), size = 0.3) +
	geom_point(mapping = aes(color = type, shape = cultivar), size = 1.5) +
	facet_wrap(~size_class, ncol = 2,
		   labeller = labeller(size_class = 
				       c("[50-100[" = "[50-100 bp[ insertions",
					 "[100-1000[" = "[100-1,000 bp[ insertions",
					 "[1000-10000[" = "[1,000-10,000 bp[ insertions",
					 "[10000+[" = "[10,000+ bp[ insertions"))) +
	x_axis +
	y_axis +
	scale_shape_discrete(name = "") +
	guides(color = FALSE,
	       shape = guide_legend(override.aes = list(size = 3))) +
	common_theme +
	theme()

# Saving as a png file
png("figure_s27.png", width = 12, height = 6, units = "in", res = 500)

grid.newpage()
# Locating the subplots in the figure, leaving some space for the "A" and "B" plot labels
pushViewport(viewport(layout = grid.layout(1, 2)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(panelA, vp = viewport(x = 0.03, width = 0.97, just = "left"))
grid.text("A", x = 0.04, y = 0.90, gp = gpar(fontsize = 24))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(panelB, vp = viewport(x = 0.03, width = 0.97, just = "left"))
grid.text("B", x = 0.04, y = 0.90, gp = gpar(fontsize = 24))
popViewport()

dev.off()

