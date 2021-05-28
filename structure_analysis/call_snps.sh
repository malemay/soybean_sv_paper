#!/bin/bash

# WARNING : some programs may have to be made available in $PATH for Platypus to work properly
# module load python/2.7
# module load htslib/1.8
# module load platypus/0.8.1.1
# Getting the path to the Platypus executable from the command line
platypus=$1

# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Necessary to avoid premature termination (allows more open files)
ulimit -S -n 4096

# DEPENDENCY : structure_analysis/platypus102_Gmax_v4_bamfiles.txt
# DEPENDENCY : Illumina reads aligned with BWA and sorted/indexed with samtools
$platypus callVariants --bamFiles=platypus102_Gmax_v4_bamfiles.txt --nCPU=8 --refFile=${refgenome} \
	--skipRegionsFile=chloroplast,mitochondrion --logFileName=platypus102_Gmax_v4.log --minMapQual=20 \
       	--minBaseQual=20 --maxVariants=10 --filterReadsWithUnmappedMates=0 --filterReadsWithDistantMates=0 \
       	--filterReadPairsWithSmallInserts=0 --output=platypus_snps.vcf

