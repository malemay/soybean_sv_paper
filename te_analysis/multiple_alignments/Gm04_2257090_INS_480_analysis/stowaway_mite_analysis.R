#!/prg/R/4.0/bin/Rscript

# Loading the required packages
library(ggplot2)
library(VariantAnnotation)
library(tidyr)

# Getting the path to the executables from the command line
executables <- commandArgs(trailingOnly = TRUE)
bcftools <- executables[1]
plink <- executables[2]
bgzip <- executables[3]
tabix <- executables[4]

# Preparing the vcf files for input by compressing and indexing
# DEPENDENCY : te_analysis/query_all.vcf
system("ln -s ../../query_all.vcf")
system(paste0(bgzip, " query_all.vcf"))
system(paste0(tabix, " query_all.vcf.gz"))


# DEPENDENCY : structure_analysis/platypus_snps.vcf
system("ln -s ../../../structure_analysis/platypus_snps.vcf")
system(paste0(bgzip, " platypus_snps.vcf"))
system(paste0(tabix, " platypus_snps.vcf.gz"))

# Using bcftools to get all SNPs ± 100 kb from the insertion location
system(paste0(bcftools, " view -f PASS,. -V indels,mnps -Ov \\
	      platypus_snps.vcf.gz \\
	      Gm04:2157090-2357090 > Gm04_2157090_2357090.vcf"))

# Using plink to calculate the LD between all pairs and output the results as a matrix to "plink.ld"
system(paste0(plink, " --vcf Gm04_2157090_2357090.vcf --r2 square0 --allow-extra-chr"))

# Creating a function that generates a data.frame for plotting LD from a plink.ld file
ld_to_df <- function(filename) {

	# Reading the LD data from file
	ld_data <- read.table(filename)
	ld_data <- as.matrix(ld_data)

	# Creating variables for the number of observations and number of positions
	n_obs <- length(ld_data)
	n_pos <- nrow(ld_data)
	stopifnot(n_pos == ncol(ld_data))

	# Initializing a data.frame to put the data in for plotting with ggplot2
	ld_df <- data.frame(r2 = numeric(n_obs), x = numeric(n_obs), y = numeric(n_obs))

	# Filling the data.frame
	ld_df$r2 <- as.numeric(ld_data)
	ld_df$x  <- rep(1:nrow(ld_data), ncol(ld_data))
	ld_df$y  <- rep(1:ncol(ld_data), each = nrow(ld_data))

	ld_df
}


# Using the SNP at position 2,257,116 in the VCF file to add lines to the LD plot
vcf_file <- read.table("Gm04_2157090_2357090.vcf")
snp_pos <- which(vcf_file[[2]] == 2257116)

ggplot(ld_to_df("plink.ld"), aes(x = x, y = y, fill = -r2)) + 
	geom_tile() +
	geom_hline(yintercept = snp_pos, color = "red") +
	geom_vline(xintercept = snp_pos, color = "red")

# The boundaries of the LD block appear to be at about 320 and 475
# Let us plot this with pink lines to see if it matches
ggplot(ld_to_df("plink.ld"), aes(x = x, y = y, fill = -r2)) + 
	geom_tile() +
	geom_hline(yintercept = snp_pos, color = "red") +
	geom_vline(xintercept = snp_pos, color = "red") +
	geom_hline(yintercept = 320, color = "pink") +
	geom_vline(xintercept = 320, color = "pink") +
	geom_hline(yintercept = 475, color = "pink") +
	geom_vline(xintercept = 475, color = "pink")

# The boundaries look about right. Let us determine the physical positions of these on Gm04
vcf_file[[2]][c(320, 475)]
# [1] 2220398 2259326

# Let us use those positions to extract a new vcf file of SNPs from the platypus file
system(paste0(bcftools, " view -f PASS,. -V indels,mnps -Ov \\
	      platypus_snps.vcf.gz \\
	      Gm04:2220398-2259326 | sed s/_all\\.sort//g > Gm04_2220398_2259326.vcf"))

# Using plink to compute the LD in that window and then plot it
system(paste0(plink, " --vcf Gm04_2220398_2259326.vcf --r2 square0 --allow-extra-chr"))

# Using the SNP at position 2,257,116 in the VCF file to add lines to the LD plot
vcf_file <- read.table("Gm04_2220398_2259326.vcf")
snp_pos <- which(vcf_file[[2]] == 2257116)
nrow(vcf_file)
# [1] 156

ggplot(ld_to_df("plink.ld"), aes(x = x, y = y, fill = -r2)) + 
	geom_tile() +
	geom_hline(yintercept = snp_pos, color = "red") +
	geom_vline(xintercept = snp_pos, color = "red")

# A few of the SNPs seem out of LD with the others, but otherwise it looks fine

# Also extracting the single position corresponding to the small insertion left by the insertion of the TE
system(paste0(bcftools, " view -Ov \\
	      platypus_snps.vcf.gz \\
	      Gm04:2257092 | sed s/_all\\.sort//g > Gm04_2257092.vcf"))

# And then also extracting the position corresponding to the TE insertion itself
system(paste0(bcftools, " view -Ov \\
	      query_all.vcf.gz \\
	      Gm04:2257090 > Gm04_2257090.vcf"))

# Reading all the files as VCF objects using VariantAnnotation::readVcf
all_snps <- readVcf("Gm04_2220398_2259326.vcf")
te_ins <- readVcf("Gm04_2257090.vcf")
te_indel <- readVcf("Gm04_2257092.vcf")

# Checking if the name of the samples are identical for te_ins and te_indel
identical(colnames(geno(te_ins)[["GT"]]), colnames(geno(te_indel)[["GT"]]))

# Now cross-tabulating the genotypes at these two loci
table(geno(te_ins)[["GT"]][1, ], geno(te_indel)[["GT"]][1, ])
#      
#       ./. 0/0 0/1 1/0 1/1
#   .     2   0   0   0   0
#   ./.   1   0   0   0   0
#   0/0   5  71   1   5   8
#   1/1   9   0   0   0   0

# It appears that the indel genotype is missing whenever the insertion is present, which is not surprising
# However, when the TE insertion is absent, then the small insertion may or may not be there
# I have verified visually with the Illumina data in IGV that heterozygous calls for the indel are most
# likely homozygous for the alternate allele. Therefore, I will consider them as such here

# I will initialize a new data.frame that classifies all samples in three classes
# according to their genotype at the insertion locus : 
# - TE present
# - TE excised
# - TE absent
genotype_df <- data.frame(sample = colnames(geno(te_ins)[["GT"]]),
			  te_ins = geno(te_ins)[["GT"]][1, ],
			  te_indel = geno(te_indel)[["GT"]][1, ])

# The genotype at the insertion locus is initialized as missing
genotype_df$genotype <- NA

# It is "present" for sample for which te_ins == "1/1"
genotype_df[genotype_df$te_ins == "1/1", "genotype"] <- "present"

# It is "excised" if the genotype at the indel site is either "1/1", "0/1" or "1/0"
genotype_df[genotype_df$te_indel %in% c("1/1", "0/1", "1/0"), "genotype"] <- "excised"

# Then it is "absent" if the genotype at the indel site is "0/0"
genotype_df[genotype_df$te_indel == "0/0", "genotype"] <- "absent"

# Removing samples for which the genotype could not be determined
genotype_df <- genotype_df[!is.na(genotype_df$genotype), ]

# Creating sub-vcfs for each of these genotypes
present_samples <- genotype_df[genotype_df$genotype == "present", "sample"]
present_vcf <- geno(all_snps[, present_samples])[["GT"]]

excised_samples <- genotype_df[genotype_df$genotype == "excised", "sample"]
excised_vcf <- geno(all_snps[, excised_samples])[["GT"]]

absent_samples <- genotype_df[genotype_df$genotype == "absent", "sample"]
absent_vcf <- geno(all_snps[, absent_samples])[["GT"]]

# Creating a function to calculate alternate allele frequencies at a site, given a character vector of genotypes
alt_freq <- function(genotypes) {
	# Creating a lookup table linking the genotypes to the number of allele counts
	lookup_table <- c("./." = 0, "0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)
	stopifnot(all(genotypes %in% names(lookup_table)))

	# Computing the allele count and then allele frequency
	allele_count <- sum(lookup_table[genotypes])
	n_alleles <- 2 * (sum(genotypes != "./."))

	return(allele_count / n_alleles)
}

# Using this function to fill a data.frame with the frequencies for the different groups
freq_df <- data.frame(site = rownames(present_vcf))
freq_df$present <- apply(present_vcf, 1, alt_freq)
freq_df$excised <- apply(excised_vcf, 1, alt_freq)
freq_df$absent <- apply(absent_vcf, 1, alt_freq)
# Removing the cases where the frequency could not be computed because of missing genotypes
freq_df <- freq_df[complete.cases(freq_df), ]

# I will reformat the data.frame so the frequencies can be plotted as tiles with the three haplotypes as rows and the SNPs as columns
freq_df$site_num <- 1:nrow(freq_df)
plotting_df <- gather(freq_df, "haplotype", "alt_freq", 2:4)
plotting_df$hap_num <- ifelse(plotting_df$haplotype == "absent", 1, ifelse(plotting_df$haplotype == "present", 2, 3))

# Now the data.frame is ready for plotting
ggplot(plotting_df, aes(x = site_num, y = hap_num, fill = alt_freq)) +
	geom_tile(height = 0.5) +
	scale_fill_distiller(type = "div", palette = "RdBu") +
	theme_bw()

# I would like to add vertical lines where presence and excised haplotypes diverge
diverging_snps <- freq_df[abs(freq_df$present - freq_df$excised) > 0.3, ]

ggplot(plotting_df, aes(x = site_num, y = hap_num, fill = alt_freq)) +
	geom_tile(height = 0.5, color = "black", size = 0.5) +
	geom_vline(data = diverging_snps, mapping = aes(xintercept = site_num)) +
	scale_fill_distiller(type = "div", palette = "RdBu") +
	theme_bw()

# Calcul du nombre de sites divergents entre les différents haplotypes
sum(abs(freq_df$present - freq_df$absent) > 0.5)
# [1] 129
sum(abs(freq_df$present - freq_df$excised) > 0.5)
# [1] 3
sum(abs(freq_df$excised - freq_df$absent) > 0.5)
# [1] 126

# Saving the data.frames to files so they can be used later
save(genotype_df, file = "genotype_df.RData")
save(freq_df, file = "freq_df.RData")
save(plotting_df, file = "plotting_df.RData")
save(diverging_snps, file = "diverging_snps.RData")

