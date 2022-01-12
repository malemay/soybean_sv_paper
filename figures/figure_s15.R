#!/usr/bin/Rscript

# Figure S15 shows the proportion of true positive and false positive deletion calls for
#  different SV calling program combinations and different SV sizes

# A sample with median F1-score will be chosen so as to be representative.
#  The benchmark data on calls filtered for the number of supporting calls
#  (DP >= 2) and the number of ALT alleles within homozygous ALT calls (AC_Hom >= 4)
#  will be used.

# Loading the ggplot2 and GenomicRanges packages
library(ggplot2)
library(GenomicRanges)

# Loading the rates for all samples so we can identify a sample to use
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/sveval_ncallers_rates.RData
load("../sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/sveval_ncallers_rates.RData")

# Finding the median f-score among the overall results
deletions <- sveval_ncallers_rates$DEL
deletions <- deletions[deletions$size_class == "all", ]
deletions <- deletions[deletions$threshold == 1, ]
deletions$f1_score <- 2 * (deletions$precision * deletions$sensitivity) / (deletions$precision + deletions$sensitivity)
# We get the median row by ordering as a function of F-score and taking the 9th row
# deletions[order(deletions$f1_score), ][9, ]
#                     size_class sensitivity precision precision_shrunk  pipeline cultivar threshold  f1_score
# CAD1049_paragraph.2        all   0.5848983 0.8080104        0.7911448 paragraph  CAD1049         1 0.6785856

# So we will be using sample CAD1049. Let us load the data for that sample
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/CAD1049_paragraph.RData
load("../sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/CAD1049_paragraph.RData")

# We get the list element for which the ncallers threshold was 1
#str(sveval_output, max.level = 3)
# List of 2
#  $ qual_ths: num [1:5] 0 1 2 3 4
#  $ eval    :List of 5
#   ..$ :List of 2
#   .. ..$ eval   :Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	5 obs. of  9 variables:
#   .. ..$ regions:List of 4
#   ..$ :List of 2
#   .. ..$ eval   :Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	5 obs. of  9 variables:
#   .. ..$ regions:List of 4
#   ..$ :List of 2
#   .. ..$ eval   :Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	5 obs. of  9 variables:
#   .. ..$ regions:List of 4
#   ..$ :List of 2
#   .. ..$ eval   :Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	5 obs. of  9 variables:
#   .. ..$ regions:List of 4
#   ..$ :List of 2
#   .. ..$ eval   :Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	5 obs. of  9 variables:
#   .. ..$ regions:List of 4
# NULL
sveval_output <- sveval_output$eval[[2]]$regions
#str(sveval_output, max.level = 1)
# List of 4
#  $ INS:List of 4
#  $ DEL:List of 4
#  $ INV:List of 4
#  $ DUP:List of 4
# NULL

# From here we use a function defined in another file to extract the relevant data
# for plotting
# DEPENDENCY : scripts/format_sveval_plotting_data.R
source("../scripts/format_sveval_plotting_data.R")
plotting_data <- format_sveval_plotting_data(sveval_output, "DEL")

# Creating the plot
figure_s15 <- ggplot(plotting_data, aes(x = size / 1000, fill = truepos)) +
    facet_wrap(~ClusterIDs, ncol = 5) + geom_histogram(bins = 30) +
    scale_x_log10(name = "Deletion size (kb)",
		  breaks = c(0.1, 1, 10, 100),
		  labels = c("0.1", "1", "10", "100")) + 
    scale_y_continuous("Number of deletions") +
    scale_fill_discrete(name = "",
			labels = c("TRUE" = "True positive calls",
				   "FALSE" = "False positive calls")) +
    theme_bw() +
    theme(text = element_text(size = 14),
	  strip.text = element_text(size = 9.5),
          legend.key.height = unit(0.02, "npc"),
          legend.position = "top",
          legend.direction = "horizontal")

# Saving as a png file
# OUTPUT : figures/figure_s15.png
png("figure_s15.png", width = 10, height = 6, units = "in", res = 500)
print(figure_s15)
dev.off()

