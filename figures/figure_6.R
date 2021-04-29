#! /prg/R/4.0/bin/Rscript

# This file contains the code to generate figure 6 of the manuscript
# the figure is split over 4 panels:
# - Panel A will show the number of polymorphic LTR elements per family as compared to the results of Tian et al. (2012)
# - Panel B will show the number of polymorphic DNA transposable elements per superfamily as compared to the results of Tian et al. (2012)
# - Panel C will show the proportion of matching nucleotides between the terminal inverted repeats of various polymorphic DNA TE SVs
# - Panel D will show the allele frequencies associated with the three different alleles at the site of a polymorphic DTT insertion

# The figure will span two columns and 3 rows arranged by the grid package.
# Panels A and B will both be on the first row while panels C and D will occupy the 2nd and 3rd rows, respectively

# Loading the required libraries
library(grid)
library(ggplot2)

# Reading in the polymorphic TEs found by Tian et al. (2012)
# DEPENDENCY : tian2012_tes.txt
tian_tes <- read.table("tian2012_tes.txt", header = FALSE, stringsAsFactors = FALSE)
names(tian_tes) <- c("accession", "chrom", "pos", "family", "category", "boundary")

# Focusing on the LTR retrotrotransposons and reformatting their family names
tian_ltr <- tian_tes[tian_tes$category %in% c("LTR-RT/copia", "LTR-RT/gypsy"), ]
tian_ltr$family <- sub("RL[CG]_", "", tian_ltr$family)

# Similarly extracting and formatting my own data
polymorphic_ltr <- polymorphic_tes[polymorphic_tes$te_type %in% c("RLC", "RLG"), ]

# Checking the extent to which the family names match across both datasets
tian_family_names <- unique(tian_ltr$family)
my_family_names <- unique(polymorphic_ltr$family)

# sum(my_family_names %in% tian_family_names) / length(my_family_names)
# [1] 0.7016575
# sum(tian_family_names %in% my_family_names) / length(tian_family_names)
# [1] 0.6135266

# Looking at the unmatched names
# my_family_names[!my_family_names %in% tian_family_names]

# Those that have had names added after their ID can be easily reformatted and matched
polymorphic_ltr$family <- sub("/.*", "", polymorphic_ltr$family)
my_family_names <- unique(polymorphic_ltr$family)

# Now let us format this data into a data.frame for plotting
family_numbers <- table(polymorphic_ltr$family)
family_numbers <- data.frame(family = names(family_numbers), count = as.numeric(family_numbers))
# Adding the counts from the Tian et al. 2021 data
tian_counts <- table(tian_ltr$family)
family_numbers$tian_count <- as.numeric(tian_counts[family_numbers$family])


# We keep only families for which we have at least 10 SV matching
# and for which there was a matching fmily name in the Tian et al. data
family_numbers <- family_numbers[complete.cases(family_numbers) & family_numbers$count >= 10, ]

# Now let us plot
figure_6a <- ggplot(family_numbers, aes(x = count, y = tian_count)) +
	geom_point() +
	scale_x_log10("Number of SVs per LTR family in this study") +
	scale_y_log10("Number of occurences of LTR family in Tian et al. (2012)") +
	theme_bw()

# Now let us see what we can learn about DNA transposable elements from the ones found by Tian et al. 2021
tian_dna <- tian_tes[grepl("^DNA", tian_tes$category), ]

# Let us also extract a data.frame of DNA transposable elements from our own results
polymorphic_dna <- polymorphic_tes[polymorphic_tes$class == "II", ]

# We will create a new column to adapt our results to their classification
polymorphic_dna$tian_classes <- ifelse(!is.na(polymorphic_dna$description), polymorphic_dna$description, polymorphic_dna$superfamily)

# We will also modify the names of the category in tian_dna to match ours
te_lookup <- c("DNA/PIF" = "PIF-Harbinge",
	       "DNA/Mutator" = "Mutator",
	       "DNA/CACTA" = "CACTA",
	       "DNA/stowaway" = "MITE/stowaway",
	       "DNA/tourist" = "MITE/tourist",
	       "DNA/hAT" = "hAT",
	       "DNA/Tc1" = "Tc1-Mariner",
	       "DNA/HELITRON" = "Helitron")

tian_dna$new_category <- te_lookup[tian_dna$category]

# Extracting vectors of counts from both datasets
my_dna_counts <- table(polymorphic_dna$tian_classes)
tian_dna_counts <- table(tian_dna$new_category)
stopifnot(identical(names(my_dna_counts), names(tian_dna_counts)))
dna_counts <- data.frame(type = names(my_dna_counts),
			 count = as.numeric(my_dna_counts),
			 tian_count = as.numeric(tian_dna_counts))

# Now plotting these results
figure_6b <- ggplot(dna_counts, aes(x = count, y = tian_count)) +
	geom_point() +
	geom_text(aes(label = type), vjust = 0, hjust = 0, nudge_x = -0.05, nudge_y = -0.08, size = 2) +
	scale_x_log10("Number of SVs per DNA TE superfamily in this study") +
	scale_y_log10("Number of occurences of DNA TE superfamily in Tian et al. (2012)") +
	theme_bw() +
	theme(text = element_text(size = 10))



# Preparing panel C
# Reading the data for plotting
# DEPENDENCY : multiple_alignments/filtered_alignments/tir_similarity_all.txt
all_means <- read.table("multiple_alignments/filtered_alignments/tir_similarity_all.txt", header = TRUE, stringsAsFactors = FALSE)

ggplot(all_means, aes(x = family, y = similarity, color = superfamily)) +
	geom_jitter(size = 2, height = 0, width = 0.15) +
	scale_y_continuous(name = "Proportion of matching nucleotides in terminal repeats") +
	scale_x_discrete(name = "Name of matching TE entry in database") +
	scale_color_discrete(name = "Superfamily",
			     labels = c("DTH" = "PIF-Harbinger (DTH)",
					"DTM" = "Mutator (DTM)",
					"DTT" = "Tc1-Mariner (DTT)")) +
theme_bw() +
theme(axis.text.x = element_text(angle = 35, 
				 vjust = 1, 
				 hjust = 1,
				 size = 7),
      legend.position = "top",
      legend.direction = "horizontal")

# Reading the data for plotting
# DEPENDENCY : multiple_alignments/filtered_alignments/tir_similarity_noref.txt
noref_means <- read.table("multiple_alignments/filtered_alignments/tir_similarity_noref.txt", header = TRUE, stringsAsFactors = FALSE)

ggplot(noref_means, aes(x = family, y = similarity, color = superfamily)) +
	geom_jitter(size = 2, height = 0, width = 0.15) +
	scale_y_continuous(name = "Proportion of matching nucleotides in terminal repeats") +
	scale_x_discrete(name = "Name of matching TE entry in database") +
	scale_color_discrete(name = "Superfamily",
			     labels = c("DTH" = "PIF-Harbinger (DTH)",
					"DTM" = "Mutator (DTM)",
					"DTT" = "Tc1-Mariner (DTT)")) +
theme_bw() +
theme(axis.text.x = element_text(angle = 35, 
				 vjust = 1, 
				 hjust = 1,
				 size = 7),
      legend.position = "top",
      legend.direction = "horizontal")

# Preparing panel D of the figure
# DEPENDENCY : multiple_alignments/Gm04_2257090_INS_480_analysis/plotting_df.RData
load("multiple_alignments/Gm04_2257090_INS_480_analysis/plotting_df.RData")
# DEPENDENCY : multiple_alignments/Gm04_2257090_INS_480_analysis/diverging_snps.RData
load("multiple_alignments/Gm04_2257090_INS_480_analysis/diverging_snps.RData")

# Adding a column for the position of the SNPs to the diverging_snps data.frame
diverging_snps$snp_pos <- sapply(strsplit(diverging_snps$site, "_"),
			      function (x) x[1])

# Plotting the ALT allele frequencies within the interval
ggplot(plotting_df, aes(x = as.factor(site_num), y = as.factor(hap_num), fill = alt_freq)) +
	geom_tile(height = 0.65) + #, color = "black", size = 0.1) +
	geom_vline(data = diverging_snps, mapping = aes(xintercept = site_num), linetype = 3, size = 0.2) +
	scale_x_discrete(name = "",
			 breaks = diverging_snps$site_num,
			 labels = diverging_snps$snp_pos) +
        scale_y_discrete(name = "Allele at TE insertion site",
			 labels = c("1" = "absent",
				    "2" = "present",
				    "3" = "excised")) +
	scale_fill_distiller(name = "Alternate allele frequency", 
			     type = "div", 
			     palette = "RdBu") +
	theme_bw() + 
	theme(legend.direction = "horizontal",
	      legend.position = "top",
	      panel.grid = element_blank(),
	      panel.border = element_blank())

