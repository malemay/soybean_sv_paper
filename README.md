# Analysis code for the paper 

## Overview

This repository contains all the code needed to reproduce the analyses presented in the paper titled "Combined use of Oxford Nanopore and Illumina sequencing yields insights into soybean structural variation biology"

As a disclaimer, readers should be aware that most of the code was reorganized and integrated into the Makefile only after analyses were performed.
Therefore, those trying to run the analyses might run into issues related to paths or software version incompatibilities.
We encourage users who encounter issues while trying to run this code to open an issue or contact the repo maintainer directly.
We believe that the code in this repository and the Makefile should still be useful to help those interested in understanding the analyses that were performed.

## Software dependencies

The following software should be installed to reproduce the analyses.
Versions used for this work are indicated in parentheses.
The path to each of the executables should be modified in the Makefile for the code to run properly.

* [AGE](https://github.com/abyzovlab/AGE) (commit 6fa60999f573998a95b6ef751b454e6719b1849d)
* [AsmVar](https://github.com/bioinformatics-centre/AsmVar) (830b1e35157d2a4151dd8462354fb7e0ab81aa0f)
* [BayesTyper](https://github.com/bioinformatics-centre/BayesTyper) (v. 1.5)
* [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) (v. 38.25)
* [bamaddrg](https://github.com/ekg/bamaddrg) (commit 3fccbf057eef21f6304fade6c306c5bb64158865)
* [bwa](https://github.com/lh3/bwa) (v. 0.7.17)
* [bcftools](https://github.com/samtools/bcftools) (commit 7cd83b71405c993d2f81fb91126abe3e44344398)
* [blastn](http://blast.ncbi.nlm.nih.gov) (v. 2.11.0+)
* [Circos](http://circos.ca) (v. 0.69-8)
* [fastStructure](https://github.com/rajanil/fastStructure) (v. 1.0)
* [FLASH](https://sourceforge.net/projects/flashpage/) (v. 1.2.11)
* [GenericRepeatFinder](https://github.com/bioinfolabmu/GenericRepeatFinder) (commit 35b1c4d6b3f6182df02315b98851cd2a30bd6201)
* [htslib](https://github.com/samtools/htslib) (v. 1.10.2)
* [KMC](https://github.com/refresh-bio/KMC) (v. 3.0.0)
* [LAST](https://gitlab.com/mcfrith/last) (v. 1047)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) (v. 7.475)
* [manta](https://github.com/Illumina/manta) (v. 1.6.0)
* [minimap2](https://github.com/lh3/minimap2) (commit c9874e2dc50e32bbff4ded01cf5ec0e9be0a53dd)
* [NanoPlot](https://github.com/wdecoster/NanoPlot) (v. 1.18.2)
* [ngmlr](https://github.com/philres/ngmlr) (v. 0.2.7)
* [paragraph](https://github.com/Illumina/paragraph) (v. 2.4a)
* [Platypus](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data) (v. 0.8.1.1)
* [PLINK](https://www.cog-genomics.org/plink2) (v. 1.90b5.3)
* [Porechop](https://github.com/rrwick/Porechop) (commit 109e437280436d1ec27e5a5b7a34ffb752176390)
* [`R` programming language](https://cran.r-project.org/) (version 4.0.3 for analyses, 3.5.0 for figures)
* [samtools](https://github.com/samtools/samtools) (commit 26d7c73c690d298c3d4f6979224933d2e2d102cf)
* [smoove](https://github.com/brentp/smoove) (v. 0.2.4)
* [SOAPdenovo2](https://github.com/aquaskyline/SOAPdenovo2) (v. 2.04)
* [Sniffles](https://github.com/fritzsedlazeck/Sniffles) (v. 1.0.11)
* [SvABA](https://github.com/walaj/svaba) (commit e2c11a9484c9fd7761ebf0bed6425a6db37591fc)
* [SVanalyzer](https://github.com/nhansen/SVanalyzer) (commit 6a18fa3d2fcc371fb1d92955b7f45c923f3b380b)
* [vcftools](https://github.com/vcftools/vcftools) (v. 0.1.16)
* [vg](https://github.com/vgteam/vg) (v. 1.23.0 "Lavello")
* [wtdbg2](https://github.com/ruanjue/wtdbg2) (commit 79334f4d92084f5f5ff81b48f6b993ae14b1d88d)

Some programs needed for reproducing analyses were modified from existing software:

* We forked `R` package [sveval](https://github.com/jmonlong/sveval) and slightly modified it to add support for benchmarking duplications and for extracting more exhaustive output.
Thsi version can be installed from [our fork](https://github.com/malemay/sveval.git) by installing the commit 65f2781cad9c1e0979c93efac41f8157a436703f on branch soybean-nanopore-svs.

* The script [scripts/addMissingPaddingGmax4.py](https://github.com/malemay/soybean_sv_paper/blob/master/scripts/addMissingPaddingGmax4.py) was adapted from [addMissingPaddingHg38.py](https://github.com/vgteam/sv-genotyping-paper/blob/master/human/misc-scripts/addMissingPaddingHg38.py) to use the soybean reference genome instead of human reference genome.
The original MIT copyright notice is included in our modified file.

## Data availability

### Sequencing data

### Reference data

### Key VCF files

### External data

## Querying the Makefile

## Citation


