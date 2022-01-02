#!/prg/R/4.0/bin/Rscript

# Figure 4 shows the results of the fastStructure analysis on SNVs and SVs

# Loading the required packages
library(ggplot2)
library(grid)

# Loading the structure data from file
load("~/analyse_nanopore/manuscript_revision/population_analysis/structure/structure_data.RData")

# A plotting function that generates a structure plot from structure results using grid
plot_structure <- function(x, scolors, vp = NULL, differences = NULL, diffcol = NULL) {

	if(!is.null(vp)) {
		pushViewport(vp)
		on.exit(upViewport(), add = TRUE)
	}

	pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = nrow(x))))

	for(i in 1:nrow(x)) {
		pushViewport(viewport(layout.pos.col = i))
		grid.rect(y = cumsum(x[i, ]), height = x[i, ], just = c("centre", "top"),
			  gp = gpar(fill = scolors, lex = 0.1))

		if(!is.null(differences) && differences[i]) {
			grid.lines(x = unit(c(0.5, 0.5), "npc"),
			   y = unit(1, "npc") + unit(c(0.1, 0.9), "lines"),
				   gp = gpar(lty = "11", lwd = 0.8))
		}

		popViewport()
	}

	grid.yaxis(gp = gpar(fontsize = 8))
	grid.text("Ancestry", x = unit(-2.5, "lines"), rot = 90, gp = gpar(fontsize = 10))

	popViewport()

	invisible(NULL)
}

# Saving the figure as Figure 4
png("figure_4.png", width = 15, height = 12, units = "cm", res = 800)

# Plotting the data for both SNPs and SVs in the same plot using grid
grid.newpage()

# Arranging a two-row layout with a plot viewport for each Structure plot
pushViewport(viewport(layout = grid.layout(ncol = 1, nrow = 2)))

snp_viewport <- plotViewport(c(0.5, 3.1 , 0.5, 0.5), name = "snps")
sv_viewport <- plotViewport(c(0.5, 3.1, 0.5, 0.5), name = "snps")

pushViewport(viewport(layout.pos.row = 1))
grid.text("A", x = 0.025, y = 0.94, gp = gpar(fontsize = 18))
plot_structure(snp_sorted, scolors = gg_color_hue(5), vp = snp_viewport)
upViewport()

pushViewport(viewport(layout.pos.row = 2))
grid.text("B", x = 0.025, y = 0.94, gp = gpar(fontsize = 18))
plot_structure(sv_sorted, scolors = gg_color_hue(5), vp = sv_viewport, differences = differences, diffcol = diff_colors)
upViewport()
dev.off()

