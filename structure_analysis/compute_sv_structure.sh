#!/bin/bash

# WARNING : python 2.7 may need to be added to $PATH for fastStructure to work
# Loading the required modules
# module load python/2.7

# Getting the path to the fastStructure executable from the command line
faststructure=$1

# Running fastStructure with K=5 as in Torkamaneh et al. 2018
# DEPENDENCY : filtered SVs prepared for input to PLINK for PCA
$faststructure -K 5 --input=sv_pca/sv_pca_input --output=sv_structure --full 

