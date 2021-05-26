#!/prg/R/4.0/bin/Rscript

# Loading the required libraries
library(ggplot2)
library(grid)

# Loading the plotting data for deletions
# DEPENDENCY : breakpoint_refinement_analysis/deletions.RData
load("../breakpoint_refinement_analysis/deletions.RData")
# DEPENDENCY : breakpoint_refinement_analysis/insertions.RData
load("../breakpoint_refinement_analysis/insertions.RData")

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
	ggplot(deletions[deletions$pipeline == "bayestyper", ], aes(x = sensitivity, y = precision)) +
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
	ggplot(insertions[insertions$pipeline == "bayestyper", ], aes(x = sensitivity, y = precision)) +
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
	#         scale_color_discrete(guide = guide_legend(override.aes = list(color = "white"))) +
	#         scale_shape_discrete(guide = guide_legend(override.aes = list(color = "white"))) +
	common_theme +
	theme()
	#         theme(legend.text = element_text(color = "white"),
	#               legend.title = element_text(color = "white"),
	#               legend.key = element_rect(fill = "white"))

# Saving as a png file
# OUTPUT : figures/figure_s17.png
png("figure_s17.png", width = 12, height = 6, units = "in", res = 500)

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

