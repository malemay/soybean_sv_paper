#!/usr/bin/Rscript

# Figure S16 shows the proportion of true positive and false positive insertion calls for
#  different SV calling program combinations and different SV sizes

# Loading the ggplot2 and GenomicRanges packages
library(ggplot2)
library(GenomicRanges)

# Loading the rates for all samples so we can identify a sample to use
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/sveval_ncallers_rates.RData
load("../sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/sveval_ncallers_rates.RData")

# Finding the median f-score among the overall results
insertions <- sveval_ncallers_rates$INS
insertions <- insertions[insertions$size_class == "all", ]
insertions <- insertions[insertions$threshold == 1, ]
insertions$f1_score <- 2 * (insertions$precision * insertions$sensitivity) / (insertions$precision + insertions$sensitivity)
# We get the median row by ordering as a function of F-score and taking the 9th row
# insertions[order(insertions$f1_score), ]
#                     size_class sensitivity precision precision_shrunk  pipeline cultivar threshold  f1_score
# CAD1070_paragraph.2        all   0.1657478 0.7171344        0.7087452 paragraph  CAD1070         1 0.2692623
# CAD1089_paragraph.2        all   0.1801718 0.7561102        0.7439956 paragraph  CAD1089         1 0.2910016
# CAD1096_paragraph.2        all   0.2522856 0.7235529        0.7107553 paragraph  CAD1096         1 0.3741234
# CAD1049_paragraph.2        all   0.2467434 0.7744807        0.7624072 paragraph  CAD1049         1 0.3742528
# CAD1002_paragraph.2        all   0.2493558 0.7561380        0.7433659 paragraph  CAD1002         1 0.3750345
# CAD1092_paragraph.2        all   0.2460361 0.7911141        0.7788698 paragraph  CAD1092         1 0.3753413
# CAD1077_paragraph.2        all   0.2544832 0.7330076        0.7220167 paragraph  CAD1077         1 0.3778022
# CAD1022_paragraph.2        all   0.2658194 0.6948877        0.6795486 paragraph  CAD1022         1 0.3845389
# CAD1015_paragraph.2        all   0.2605210 0.7483801        0.7361268 paragraph  CAD1015         1 0.3864973
# CAD1056_paragraph.2        all   0.2680957 0.7116722        0.6985463 paragraph  CAD1056         1 0.3894724
# CAD1087_paragraph.2        all   0.2707196 0.7104459        0.6971246 paragraph  CAD1087         1 0.3920473
# CAD1018_paragraph.2        all   0.2611572 0.8089443        0.7983949 paragraph  CAD1018         1 0.3948441
# CAD1074_paragraph.2        all   0.2769189 0.7358364        0.7233146 paragraph  CAD1074         1 0.4024013
# CAD1065_paragraph.2        all   0.2831057 0.7143872        0.7034305 paragraph  CAD1065         1 0.4055108
# CAD1064_paragraph.2        all   0.3011564 0.7048236        0.6957328 paragraph  CAD1064         1 0.4220008
# CAD1052_paragraph.2        all   0.3000698 0.7725278        0.7613314 paragraph  CAD1052         1 0.4322446
# CAD1010_paragraph.2        all   0.3222873 0.7300100        0.7215321 paragraph  CAD1010         1 0.4471606

# Would probably be best to choose the median again. This can be easily justified in the paper,
#  especially since this is for supplementary material anyway
# We get the median row by ordering as a function of F-score and taking the 9th row
# insertions[order(insertions$f1_score), ][9, ]
#                     size_class sensitivity precision precision_shrunk  pipeline cultivar threshold  f1_score
# CAD1015_paragraph.2        all    0.260521 0.7483801        0.7361268 paragraph  CAD1015         1 0.3864973

# So we will be using sample CAD1015. Let us load the data for that sample
# DEPENDENCY : sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/CAD1015_paragraph.RData
load("../sv_genotyping/illumina_svs/sveval_benchmarks/ncallers_RData/CAD1015_paragraph.RData")

# We get the list element for which the ncallers threshold was 1
sveval_output <- sveval_output$eval[[2]]$regions

# From here we use a function defined in another file to extract the relevant data
# for plotting
# DEPENDENCY : scripts/format_sveval_plotting_data.R
source("../scripts/format_sveval_plotting_data.R")
plotting_data <- format_sveval_plotting_data(sveval_output, "INS")

# Creating the plot
figure_s16 <- ggplot(plotting_data, aes(x = size, fill = truepos)) +
    facet_wrap(~ClusterIDs) + geom_histogram(bins = 30) +
    scale_x_log10(name = "Insertion size (bp)") + 
    scale_y_continuous("Number of insertions") +
    scale_fill_discrete(name = "",
			labels = c("TRUE" = "True positive calls",
				   "FALSE" = "False positive calls")) +
    theme_bw() +
    theme(text = element_text(size = 14),
          legend.key.height = unit(0.02, "npc"),
          legend.position = "top",
          legend.direction = "horizontal")

# Also saving as a png file
# OUTPUT : figures/figure_s16.png
png("figure_s16.png", width = 6, height = 6, units = "in", res = 500)
print(figure_s16)
dev.off()

