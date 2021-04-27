#!/prg/R/4.0/bin/Rscript

# Figure 4 of the manuscript

# This figure shows the 4 first principal components obtained for
# both the PCA conducted using SNPs and using SVs
#
# Points will be colored according to their population assignment
# with fastStructure; points with less than 0.5 as their maximum
# q-value will be shown with a different symbol to indicate admixture

# Loading the required packages
library(ggplot2)
library(grid)

# Creating variables for the directories containing the analyses based on SVs and SNPs
sv_dir <- "/home/malem420/WGS_data/bbduk_trimmed/bwa_alignment_Gmax_v4/paragraph_all/paragraph_illumina/pca_analysis/pca_final"
snp_dir <- "/home/malem420/WGS_data/bbduk_trimmed/structure/Gmax_v4"

# First checking that the samples were in the same order for both analyses
# DEPENDENCY : snp_samples.txt, CAD.ped
system(paste0("cut -f1 ", snp_dir, "/CAD.ped > snp_samples.txt"))
# DEPENDENCY : sv_samples.txt, paragraph.ped
system(paste0("cut -f1 ", sv_dir, "/paragraph.ped > sv_samples.txt"))
snp_samples <- scan("snp_samples.txt", what = character(), quiet = TRUE, sep = "\n")
sv_samples <- scan("sv_samples.txt", what = character(), quiet = TRUE, sep = "\n")
stopifnot(identical(substr(snp_samples, 1, 7), sv_samples))

# Reading the eigenvalues for the SV analysis and computing the percentage explained variance
# DEPENDENCY : paragraph_pca.eigenval
sv_eigenval <- read.table(paste0(sv_dir, "/paragraph_pca.eigenval"), header = FALSE)
sv_explained <- sv_eigenval[[1]]/sum(sv_eigenval[[1]])

# Doing the same thing for the SNP analysis
# DEPENDENCY : CAD_pca.eigenval
snp_eigenval <- read.table(paste0(snp_dir, "/CAD_pca.eigenval"), header = FALSE)
snp_explained <- snp_eigenval[[1]]/sum(snp_eigenval[[1]])

# Reading the principal components for the SV analysis
# DEPENDENCY : paragraph_pca.eigenvec
sv_pca <- read.table(paste0(sv_dir, "/paragraph_pca.eigenvec"), header = FALSE, stringsAsFactors = FALSE)
names(sv_pca) <- c("ind", "ind2", paste0("PC", as.character(1:20)))

# Doing the same thing for the SNP analysis
# DEPENDENCY : CAD_pca.eigenvec
snp_pca <- read.table(paste0(snp_dir, "/CAD_pca.eigenvec"), header = FALSE, stringsAsFactors = FALSE)
names(snp_pca) <- c("ind", "ind2", paste0("PC", as.character(1:20)))

# Just another sanity check
stopifnot(identical(sv_pca[[1]], substr(snp_pca[[1]], 1, 7)))

# Reading the fastStructure Q matrix obtained from the SNP analysis
# DEPENDENCY : CADtest.5.meanQ
qmat <- read.table("/home/malem420/WGS_data/bbduk_trimmed/structure/Gmax_v4/CADtest.5.meanQ", header = FALSE)

# Adding a column for the population assignment in the pca data.frames
sv_pca$pop <- as.factor(apply(qmat, 1, which.max))
snp_pca$pop <- as.factor(apply(qmat, 1, which.max))

# Adding a column indicating whether the maximum Q-value for a given sample is over 0.6
sv_pca$admixed <- apply(qmat, 1, max) < 0.6
snp_pca$admixed <- apply(qmat, 1, max) < 0.6

# Creating a common theme for the plots
common_theme <- theme_bw()

# A function that generates a plot for a given PCA data.frame and principal components
pca_plot <- function(data, pc1, pc2) {
	ggplot(data[!data$admixed, ], mapping = aes_string(x = pc1, y = pc2, color = "pop")) +
		geom_point() +
		geom_point(data[data$admixed, ], mapping = aes_string(x = pc1, y = pc2, color = "pop"), shape = 8)
}

# Creating plots for the first 4 principal components for the SV dataset
sv_pc1_pc2 <- 
	pca_plot(sv_pca, "PC1", "PC2") + 
	scale_x_continuous(name = "") +
	scale_y_continuous(name = paste0("PC2 (", sprintf("%.1f", round(sv_explained[2], 3) * 100), " %)")) +
	guides(color = FALSE) +
	common_theme

sv_pc1_pc3 <- 
	pca_plot(sv_pca, "PC1", "PC3") + 
	scale_x_continuous(name = "") +
	scale_y_continuous(name = paste0("PC3 (", sprintf("%.1f", round(sv_explained[3], 3) * 100), " %)")) +
	guides(color = FALSE) +
	common_theme

sv_pc1_pc4 <- 
	pca_plot(sv_pca, "PC1", "PC4") + 
	scale_x_continuous(name = "") +
	scale_y_continuous(name = paste0("PC4 (", sprintf("%.1f", round(sv_explained[4], 3) * 100), " %)")) +
	guides(color = FALSE) +
	common_theme

# Creating the plots for the first 4 principal components for the SNP dataset
snp_pc1_pc2 <- 
	pca_plot(snp_pca, "PC1", "PC2") + 
	scale_x_continuous(name = "") +
	scale_y_continuous(name = paste0("PC2 (", sprintf("%.1f", round(snp_explained[2], 3) * 100), " %)")) +
	guides(color = FALSE) +
	common_theme

snp_pc1_pc3 <- 
	pca_plot(snp_pca, "PC1", "PC3") + 
	scale_x_continuous(name = "") +
	scale_y_continuous(name = paste0("PC3 (", sprintf("%.1f", round(snp_explained[3], 3) * 100), " %)")) +
	guides(color = FALSE) +
	common_theme

snp_pc1_pc4 <- 
	pca_plot(snp_pca, "PC1", "PC4") + 
	scale_x_continuous(name = "") +
	scale_y_continuous(name = paste0("PC4 (", sprintf("%.1f", round(snp_explained[4], 3) * 100), " %)")) +
	guides(color = FALSE) +
	common_theme

# Arranging all the plots together in a single image for publication

# tiff("/home/malem420/sv_manuscript/figures/figure_4.tiff", width = 6, height = 9, units = "in", res = 300, compression = "lzw")
# grid.newpage()
# 
# # Pushing a first viewport with a layout such that SVs will be to the left, and SNPs to the right
# pushViewport(viewport(layout = grid.layout(1, 2)))
# pushViewport(viewport(layout.pos.col = 1))
# # And now a second viewport that will be used to print the SV plots, but leaves some space for margin annotations
# pushViewport(viewport(x = 0.08, y = 0.02, width = 0.92, height = 0.98, just = c("left", "bottom")))
# # The next one splits the space for each of the three PC plots
# pushViewport(viewport(layout = grid.layout(3, 1)))
# print(sv_pc1_pc2, vp = viewport(layout.pos.row = 1))
# print(sv_pc1_pc3, vp = viewport(layout.pos.row = 2))
# print(sv_pc1_pc4, vp = viewport(layout.pos.row = 3))
# # Going up a viewport and printing the x-axis name and the panel label
# popViewport()
# grid.text(paste0("PC1 (", sprintf("%.1f", round(sv_explained[1], 3) * 100), " %)"), 
# 	  x = 0.6, y = 0)
# # Going up yet another viewport
# popViewport()
# grid.text("A", x = 0.07, y = 0.97, gp = gpar(fontsize = 24))
# popViewport()
# 
# # Now doing the sme thing for the SNP PCA
# pushViewport(viewport(layout.pos.col = 2))
# # And now a second viewport that will be used to print the SNP plots, but leaves some space for margin annotations
# pushViewport(viewport(x = 0.08, y = 0.02, width = 0.92, height = 0.98, just = c("left", "bottom")))
# # The next one splits the space for each of the three PC plots
# pushViewport(viewport(layout = grid.layout(3, 1)))
# print(snp_pc1_pc2, vp = viewport(layout.pos.row = 1))
# print(snp_pc1_pc3, vp = viewport(layout.pos.row = 2))
# print(snp_pc1_pc4, vp = viewport(layout.pos.row = 3))
# # Going up a viewport and printing the x-axis name and the panel label
# popViewport()
# grid.text(paste0("PC1 (", sprintf("%.1f", round(snp_explained[1], 3) * 100), " %)"), 
# 	  x = 0.6, y = 0)
# # Going up yet another viewport
# popViewport()
# grid.text("B", x = 0.07, y = 0.97, gp = gpar(fontsize = 24))
# popViewport()
# 
# dev.off()


# Now we create another version with just the first two principal components 
png("figure_4.png", width = 6, height = 6, units = "in", res = 500)

grid.newpage()

# Pushing a first viewport with a layout such that SVs will be to the left, and SNPs to the right
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1))
# And now a second viewport that will be used to print the SV plots, but leaves some space for margin annotations
pushViewport(viewport(x = 0.04, y = 0.02, width = 0.92, height = 0.98, just = c("left", "bottom")))
# The next one splits the space for each of the three PC plots
pushViewport(viewport(layout = grid.layout(2, 1)))
print(sv_pc1_pc2, vp = viewport(layout.pos.row = 1))
print(sv_pc1_pc3, vp = viewport(layout.pos.row = 2))
# Going up a viewport and printing the x-axis name and the panel label
popViewport()
grid.text(paste0("PC1 (", sprintf("%.1f", round(sv_explained[1], 3) * 100), " %)"), 
	  x = 0.6, y = 0)
# Going up yet another viewport
popViewport()
grid.text("A", x = 0.06, y = 0.97, gp = gpar(fontsize = 20))
popViewport()

# Now doing the sme thing for the SNP PCA
pushViewport(viewport(layout.pos.col = 2))
# And now a second viewport that will be used to print the SNP plots, but leaves some space for margin annotations
pushViewport(viewport(x = 0.08, y = 0.02, width = 0.92, height = 0.98, just = c("left", "bottom")))
# The next one splits the space for each of the three PC plots
pushViewport(viewport(layout = grid.layout(2, 1)))
print(snp_pc1_pc2, vp = viewport(layout.pos.row = 1))
print(snp_pc1_pc3, vp = viewport(layout.pos.row = 2))
# Going up a viewport and printing the x-axis name and the panel label
popViewport()
grid.text(paste0("PC1 (", sprintf("%.1f", round(snp_explained[1], 3) * 100), " %)"), 
	  x = 0.6, y = 0)
# Going up yet another viewport
popViewport()
grid.text("B", x = 0.10, y = 0.97, gp = gpar(fontsize = 20))
popViewport()

dev.off()

