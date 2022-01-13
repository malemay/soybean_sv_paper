#!/usr/bin/Rscript

# Figure s2 comprises two panels:
#  - Panel A shows the distribution of heterozygosity rates as a historgram
#  - Panel B shows the relationship between heterozygosity and sequencing depth

# Heterozygosity rates are based on the filtered set of SNVs, and sequencing depth on mapped Illumina data

# A plot of heterozygosity rates based on SNPs in the Canadian population
# DEPENDENCY : structure_analysis/snp_heterozygosity_rates.txt
hetrates <- read.table("../structure_analysis/snp_heterozygosity_rates.txt", header = TRUE, stringsAsFactors = FALSE)

hetrates$hetrate <- 1 - (hetrates$O.HOM. / hetrates$N_SITES)
hetrates$INDV <- substr(hetrates$INDV, 1, 7)

# DEPENDENCY : tables/table_s7.csv
sequencing_depth <- read.table("../tables/table_s7.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)

stopifnot(identical(hetrates$INDV, sequencing_depth$ID))

png("figure_s2.png", width = 480, height = 960, res = 200, pointsize = 8)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
hist(hetrates$hetrate, xlab = "Heterozygosity rate", ylab = "Number of cultivars", breaks = 25, main = "", cex.axis = 0.9)
mtext("A", side = 3, at = -0.135, cex = 2)
plot(sequencing_depth$depth, hetrates$hetrate, xlab = "Sequencing depth (X)", ylab = "Heterozygosity rate",
     axes = FALSE, xlim = c(0, 25), ylim = c(0, 0.4), cex = 0.9, cex.axis = 0.8)
box()
axis(side = 1, at = seq(0, 25, 5), cex.axis = 0.9)
axis(side = 2, at = seq(0, 0.4, 0.05), cex.axis = 0.9)
mtext("B", side = 3, at = -9, cex = 2)
dev.off()

