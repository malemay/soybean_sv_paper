#!/prg/R/4.0/bin/Rscript

# Code to create Figure S4 of the manuscript
# This figure is a scatter plot of the precision of dupliction genotyping
#  as a function of Oxford Nanopore sequencing depth
# The ojective is to show that samples that were sequenced using Nanopore
#  have higher precision of duplication genotyping

# Loading the ggplot2 and grid packages
library(ggplot2)
library(grid)

# Loading the sequencing depth data
# DEPENDENCY : average_depth.RData
load("/home/malem420/sv_manuscript/depth_distributions/average_depth.RData")

# Loading the genotyping precision and sensitivity rates
# DEPENDENCY : sveval_nogeno_rates.RData
load("/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_genotyping/svmerged_variants/paragraph/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData")

# We also load the correspondence between line names and their CAD IDs
# DEPENDENCY : line_ids.txt
line_ids <- read.table("/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_genotyping/svmerged_variants/paragraph/sveval_benchmarks/line_ids.txt", 
		       header = TRUE, stringsAsFactors = FALSE)

# Now we can extract and link the data through the IDs
line_ids$sequencing_depth <- average_depth[match(line_ids$name, average_depth$sample), "depth"]

# Checking whether there is a link between Oxford Nanopore sequencing depth and precision all SV types
# First writing a function that will extract the data from the rates object
get_precision <- function(line_ids, precision_data, svtype) {
	precision_data <- precision_data[[svtype]]
	precision_data <- precision_data[precision_data$size_class == "all" & precision_data$threshold == 2, ]
	precision <- precision_data$precision
	names(precision) <- precision_data$cultivar

	line_ids$precision <- precision[line_ids$id]
	return(line_ids)
}

# Creating a common theme for all panels
common_theme <- theme_bw() + 
	theme(panel.grid.minor = element_blank(),
	      text = element_text(size = 6))

# Creating the plot object for deletions
del_plot <- ggplot(get_precision(line_ids, sveval_nogeno_rates, "DEL"), aes(x = sequencing_depth, y = precision)) +
	geom_point(size = 1) +
	scale_x_continuous(name = "Average sequencing depth (X)") +
	scale_y_continuous(name = "Precision",
			   limits = c(0.70, 0.87)) +
	ggtitle("Deletions") +
	common_theme

# Creating the plot object for insertions
ins_plot <- ggplot(get_precision(line_ids, sveval_nogeno_rates, "INS"), aes(x = sequencing_depth, y = precision)) +
	geom_point(size = 1) +
	scale_x_continuous(name = "Average sequencing depth (X)") +
	scale_y_continuous(name = "Precision") +
	ggtitle("Insertions") +
	common_theme

# Creating the plot object for duplications
dup_plot <- ggplot(get_precision(line_ids, sveval_nogeno_rates, "DUP"), aes(x = sequencing_depth, y = precision)) +
	geom_point(size = 1) +
	scale_x_continuous(name = "Average sequencing depth (X)") +
	scale_y_continuous(name = "Precision",
			   breaks = seq(0.10, 0.22, 0.04)) +
	ggtitle("Duplications") +
	common_theme

# Creating the plot object for inversions
inv_plot <- ggplot(get_precision(line_ids, sveval_nogeno_rates, "INV"), aes(x = sequencing_depth, y = precision)) +
	geom_point(size = 1) +
	scale_x_continuous(name = "Average sequencing depth (X)") +
	scale_y_continuous(name = "Precision") +
	ggtitle("Inversions") +
	common_theme

# Saving to disk as a png file
png("figure_s4.png", width = 3, height = 3, units = "in", res = 500)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
del_vp <- viewport(layout.pos.col = 1, layout.pos.row = 1)
ins_vp <- viewport(layout.pos.col = 2, layout.pos.row = 1)
inv_vp <- viewport(layout.pos.col = 1, layout.pos.row = 2)
dup_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2)
print(del_plot, vp = del_vp)
print(ins_plot, vp = ins_vp)
print(inv_plot, vp = inv_vp)
print(dup_plot, vp = dup_vp)
grid.text("A", x = 0.05, y = 0.95, vp = del_vp)
grid.text("B", x = 0.05, y = 0.95, vp = ins_vp)
grid.text("C", x = 0.05, y = 0.95, vp = dup_vp)
grid.text("D", x = 0.05, y = 0.95, vp = inv_vp)
dev.off()


