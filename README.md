# Code for the analysis of structural variation in soybean

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

The [breakpoint refinement pipeline](https://github.com/malemay/breakpoint_refinement) should be installed under `scripts/breakpoint_refinement`.

Some programs needed for reproducing analyses were modified from existing software:

* We forked `R` package [sveval](https://github.com/jmonlong/sveval) and slightly modified it to add support for benchmarking duplications and for extracting more exhaustive output.
This version can be installed from [our fork](https://github.com/malemay/sveval.git) by installing the commit 65f2781cad9c1e0979c93efac41f8157a436703f on branch soybean-nanopore-svs.

* The script [scripts/addMissingPaddingGmax4.py](https://github.com/malemay/soybean_sv_paper/blob/master/scripts/addMissingPaddingGmax4.py) was adapted from [addMissingPaddingHg38.py](https://github.com/vgteam/sv-genotyping-paper/blob/master/human/misc-scripts/addMissingPaddingHg38.py) to use the soybean reference genome instead of human reference genome.
The original MIT copyright notice is included in our modified file.

## Data availability

### Sequencing data

The Illumina data used in this project is available from the [SRA](https://www.ncbi.nlm.nih.gov/sra) using the study accession number [SRP094720](https://www.ncbi.nlm.nih.gov/sra/?term=SRP094720).
These data should be placed under `illumina_data/raw_fastq/` to reproduce the analyses.

the Oxford Nanopore data generated by this project is available from the [SRA](https://www.ncbi.nlm.nih.gov/sra) using the study accession number [SRP331097](https://www.ncbi.nlm.nih.gov/sra/?term=SRP331097).
These data should be placed under `nanopore_data/` to reproduce the analyses.

### Reference data

The following datasets are available from the Web and should be added to the repository to reproduce the analyses:

* The SoyTEdb fasta file (`SoyBase_TE_Fasta.txt`) can be downloaded from [SoyBase](https://www.soybase.org/soytedb/) and should be placed under `te_analysis/te_database/` to reproduce the analyses.
* The non-reference tranposable elements found by Tian et al. (2012) can be downloaded from the supplementary material to [their paper](https://doi.org/10.1105/tpc.112.103630).
The data can be converted to a text file and saved under `te_analysis/tian2012_tes.txt`.
* The reference genome sequence and annotation of soybean cultivar Williams82, assembly version 4 can be downloaded from [Phytozome](https://phytozome-next.jgi.doe.gov/).
The files needed (`Gmax_508_v4.0.fa`, `Gmax_508_Wm82.a4.v1.gene_exons.gff3`, `Gmax_508_Wm82.a4.v1.gene.gff3`, `Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3`) should be placed under `refgenome/`.
* The soybean chloroplast and mitochondrion genome sequences can be downloaded from [SoyBase](https://www.soybase.org/GlycineBlastPages/blast_descriptions.php).
These should be concatenated together and placed under `refgenome/bt_decoy_sequences.fasta`.
They should also be concatenated to `Gmax_508_v4.0.fa` and placed under the name `refgenome/Gmax_508_v4.0_mit_chlp.fasta`


### Data generated by the analysis

Several of the VCF files generated by the analysis as well as the result from the permutation test on the overlap between SVs and genic features are available on [Figshare](https://doi.org/10.6084/m9.figshare.15127730.v1).

## Querying the Makefile

We used [GNU Make](https://www.gnu.org/software/make/) to describe the dependencies among our scripts and data through a Makefile.
In theory, the Makefile should allow running all the analyses that were done in the paper in the proper order, given that all the sequencing and reference data are made available.
In practice, our Makefile is intended as a tool to query this repository to understand what scripts should be run and in what order to obtain a particular result.
Here, we give a short introduction for people who are not yet familiar with Make so they can query our Makefile effectively.

GNU make should be installed by default on many Linux distributions.
If not, please visit [their website](https://www.gnu.org/software/make/) for download and installation.

Make describes dependencies using a set of targets which depend on a list of prerequisites, and includes for each target a so-called recipe of shell commands used to create the target from the prerequisites.
To get a list of the available targets in the Makefile, simply type the following command while in the top-level directory of this repository:

	make list

Each of these targets can be given as a argument to the `make` command to launch the commands required to create the target.
For example, the following command would run all the analyses used to make the paper:

	make all

We do not recommend running this command as such given the high computing requirements.
However, the list of all commands that would be run if the command were to be launched can be obtained with the `-n` option:

	make -n all

To get all the commands lancuhed to create Figure 1 from scratch, the following command can be used:

	make -n figures/figure_1.png

Make automatically determines which commands need to be run based on the last time when each target was updated.
If you want to trick Make into thinking that all targets were properly created, you can use the `-t` option to `touch` each target and update its timestamp instead of running the commands:

	make -t all

If this command runs properly, then running `make all` should print `make: Nothing to be done for 'all'.`
If you want to for the list of commands to be printed even though the target is up to date, then you can add the `-B` option:

	make -Bn all

You can simulate what commands would actually be run from a given point in the analysis by modifying the timestamp of a code or data file of interest and running `make -n` again.
For example, the code below would list all the commands needed to update Figure 3 after the code to call SNVs using playtpus has been modified.

	touch structure_analysis/call_snps.sh
	make -n figures/figure_3.png

With these tools in hand, you should be able to effectively query the Makefile and understand the analysis workflow that we used.

## Citation

If you use part of this code for your analyses, please cite:

(to be added later)

