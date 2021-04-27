#!/prg/R/4.0/bin/Rscript

# Code for the Figure S8 of the manuscript
# This figure is divided into two parts:
# - Panel A shows the number of variants of various sizes in the range -250 to 250
#   for various tools as histograms
# - Panel B shows the distribution of deletion sizes on a logorithmic scale for
#   various tools as histograms

# Loading the required libraries
library(ggplot2)
library(grid)

# Reading in the data that will be used for plotting
# DEPENDENCY : size_distribution.tsv
sv_sizes <- read.table("/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_genotyping/svmerged_variants/size_distribution.tsv",
		       header = TRUE, stringsAsFactors = FALSE)
# nrow(sv_sizes)
# [1] 140597

# Removing the few SVs that might be shorter than 50 bp
sv_sizes <- sv_sizes[abs(sv_sizes$size) >= 50, ]
# nrow(sv_sizes)
# [1] 138000

# Creating a common theme for the two panels
common_theme <- 
	theme_bw() +
	theme(text = element_text(size = 14),
	      axis.text = element_text(size = 10),
	      panel.grid.minor = element_blank(),
	      strip.text = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")))

# Creating panel A
panelA <- ggplot(sv_sizes[abs(sv_sizes$size) <= 250 & sv_sizes$svtype %in% c("INS", "DEL"), ], aes(x = size)) +
	geom_histogram(binwidth = 10, fill = "skyblue", color = "black", size = 0.15) +
	facet_wrap(~program, ncol = 1) +
	scale_x_continuous(name = "SV size (bp)") +
	ylab("Number of SVs") +
	common_theme

# Creating panel B
panelB <- ggplot(sv_sizes[sv_sizes$svtype == "DEL", ], aes(x = abs(size))) +
	geom_histogram(fill = "indianred", color = "black", bins = 50, size = 0.15) +
	facet_wrap(~program, ncol = 1) +
	scale_x_log10(name = "Deletion size (bp)",
		      breaks = c(100, 10000, 1000000)) +
	ylab("Number of deletions") +
	common_theme +
	theme(panel.grid.minor = element_line(),
	      panel.grid.minor.y = element_blank())

# Assembling both panels in a single figure and saving to "figure_s8.png"
png("figure_s8.png", width = 6, height = 6, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(panelA, vp = viewport(x = 0.05, width = 0.95, just = "left"))
grid.text("A", x = 0.07, y = 0.97, gp = gpar(fontsize = 24))
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(panelB, vp = viewport(x = 0.05, width = 0.95, just = "left"))
grid.text("B", x = 0.07, y = 0.97, gp = gpar(fontsize = 24))
dev.off()

