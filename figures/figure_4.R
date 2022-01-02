#!/prg/R/4.0/bin/Rscript

# Figure 4 shows the results of the fastStructure analysis on SNVs and SVs

# Loading the required packages
library(combinat)
library(ggplot2)
library(grid)

# Reading the results of the structure analysis from SNPs
# DEPENDENCY : structure_analysis/structure.5.meanQ
snp_structure <- read.table("../structure_analysis/structure.5.meanQ")
snp_structure <- as.matrix(snp_structure)
rowSums(snp_structure)

# Reading the results of the structure analysis from SVs
# DEPENDENCY : structure_analysis/sv_structure.5.meanQ
sv_structure <- read.table("../structure_analysis/sv_structure.5.meanQ")
sv_structure <- as.matrix(sv_structure)
rowSums(sv_structure)

# We need to find a permutation that makes the population numbers (1 to 5) equivalent for both analyses
# Writing a function that finds the permutation of the second data.frame which minimizes the difference
find_perm <- function(x, y) {
	# Checking that the dimensions of the two objects are the same
	stopifnot(identical(dim(x), dim(y)))

	# Generating a matrix with all permutations as rows
	permutations <- as.matrix(t(unname(as.data.frame(permn(1:ncol(y))))))

	# Creating a vector that will hold the sum of distances
	distances <- numeric(nrow(permutations))

	# Iterating over the possible permutations
	for(i in 1:nrow(permutations)) {
		distances[i] <- sum(abs(x - y[, permutations[i, ]]))
	}

	# Returning the y matrix with the best permutation
	y[, permutations[which.min(distances), ]]
}

sv_structure_perm <- find_perm(snp_structure, sv_structure)
stopifnot(all(rowSums(sv_structure_perm) == rowSums(sv_structure)))

sum(apply(snp_structure, 1, which.max) != apply(sv_structure_perm, 1, which.max))

# Sorting the samples within the first matrix according to their population estimates
sample_sort <- order(-apply(snp_structure, 1, which.max), snp_structure[, 1], snp_structure[, 2], snp_structure[, 3], snp_structure[, 4], snp_structure[, 5], decreasing = TRUE)
snp_sorted <- snp_structure[sample_sort, ]
sv_sorted  <- sv_structure_perm[sample_sort, ] 

# A function to imitate the color selection in ggplot2
# Taken from Stack Overflow: https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

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

# A vector indicating whether the population with the highest value is the same for both analyses
differences <- apply(snp_sorted, 1, which.max) != apply(sv_sorted, 1, which.max)
diff_colors <- ifelse(apply(snp_sorted, 1, max) < 0.6, "black", "red")

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

