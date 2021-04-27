#!/prg/R/4.0/bin/Rscript

# Code to create Figure S4 of the manuscript
# This figure is a scatter plot of the precision of dupliction genotyping
#  as a function of Oxford Nanopore sequencing depth
# The ojective is to show that samples that were sequenced using Nanopore
#  have higher precision of duplication genotyping

# Loading the ggplot2 package
library(ggplot2)

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

duplication_data <- sveval_nogeno_rates$DUP
duplication_data <- duplication_data[duplication_data$size_class == "all" & duplication_data$threshold == 2, ]
dup_precision <- duplication_data$precision
names(dup_precision) <- duplication_data$cultivar

line_ids$precision <- dup_precision[line_ids$id]

# Creating the plot object
figure_s4 <- ggplot(line_ids, aes(x = sequencing_depth, y = precision)) +
	geom_point(size = 1) +
	scale_x_continuous(name = "Average Oxford Nanopore sequencing depth (X)",
			   limits = c(5, 25),
			   breaks = seq(5, 25, 5)) +
	scale_y_continuous(name = "Precision of duplication genotyping",
			   breaks = seq(0.10, 0.22, 0.02)) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
	      text = element_text(size = 8))

# Saving to disk as a png file
png("figure_s4.png", width = 3, height = 3, units = "in", res = 500)
print(figure_s4)
dev.off()

