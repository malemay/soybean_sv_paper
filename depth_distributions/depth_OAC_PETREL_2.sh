#!/bin/bash

#SBATCH -D /home/malem420/sv_manuscript/depth_distributions
#SBATCH -J depth_OAC_PETREL_2
#SBATCH -o depth_OAC_PETREL_2_%j.out
#SBATCH -c 1
#SBATCH -p soyagen
#SBATCH -A soyagen
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marc-andre.lemay.2@ulaval.ca
#SBATCH --time=04:00:00
#SBATCH --mem=2000

# Computing the sequencing depth distribution on OAC_PETREL_2
# The sequencing depth for this sample had not been computed before because it had
# only been computed for OAC_PETREL_1 and OAC_PETREL_2 merged together, but we
# have since dropped OAC_PETREL_1

# Running the bash/awk script on the bam file
input_file=/home/malem420/analyse_nanopore/OAC_PETREL_2/guppy_4.0.11/ngmlr_alignment/OAC_PETREL_2_porechopped_aligned.sort.bam
/home/malem420/scripts/depth_distribution.sh $input_file OAC_PETREL_2_depth.txt

