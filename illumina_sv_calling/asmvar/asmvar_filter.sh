#!/bin/bash

# Creating a variable for the bcftools executable
bcftools=$1
bayestypertools=$2

# --- 1: Gathering and filtering all the vcf files for a given sample into files following pattern asmvar_filtering/CAD????_filtered.vcf

# Creating a variable to hold the names of the 20 chromosomes
chrlist=$(seq -w 1 20 | xargs -I {} echo -n "Gm{} ")

# Creating a variable for the IDs of all the cultivars
# DEPENDENCY : utilities/all_lines.txt
lines=$(cat ../../utilities/all_lines.txt)

# The i variable iterates over all cultivars
for i in $lines
do
    # The j variable iterates over all chromosomes
    for j in $chrlist
    do
	    # DEPENDENCY : VCF files output by AsmVar
	    # Creating a variable for the name of the file
	    ij_file=asmvar_calling/${i}/asmvar_results_${j}.vcf
	    # Adding a header with all contig names based on the reference, and saving the result to a new file
	    # DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
	    $bcftools reheader --fai ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai $ij_file > asmvar_filtering/${i}_${j}_header.vcf
    done
  
    # Concatenating the files for all chromosomes together, removing records without a "." FILTER tag, and removing SNPs
    $bcftools concat -Ou $(ls asmvar_filtering/${i}_Gm??_header.vcf) | \
	    $bcftools view --apply-filters .,PASS -Ou - | \
	    $bcftools view -Ov --exclude-types snps - > asmvar_filtering/${i}_filtered.vcf

    # Deleting the temporary files with the added headers
    rm asmvar_filtering/${i}_Gm??_header.vcf
done

# --- 2: Normalizing the allelic representation in the VCF files of each sample

# Looping over the filtered vcf files for all lines
for i in $lines
do
	# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
	$bcftools norm -Ov -f ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta asmvar_filtering/${i}_filtered.vcf > asmvar_filtering/${i}_norm.vcf
done

# --- 3: Combining the variants from all samples into a single vcf file

ulimit -n 2048

# Moving to the directory for simplicity
cd asmvar_filtering

$bayestypertools combine -o asmvar_combined -v CAD1001:CAD1001_norm.vcf,CAD1002:CAD1002_norm.vcf,CAD1003:CAD1003_norm.vcf,CAD1004:CAD1004_norm.vcf,CAD1005:CAD1005_norm.vcf,CAD1006:CAD1006_norm.vcf,CAD1007:CAD1007_norm.vcf,CAD1008:CAD1008_norm.vcf,CAD1009:CAD1009_norm.vcf,CAD1010:CAD1010_norm.vcf,CAD1011:CAD1011_norm.vcf,CAD1012:CAD1012_norm.vcf,CAD1013:CAD1013_norm.vcf,CAD1014:CAD1014_norm.vcf,CAD1015:CAD1015_norm.vcf,CAD1016:CAD1016_norm.vcf,CAD1017:CAD1017_norm.vcf,CAD1018:CAD1018_norm.vcf,CAD1019:CAD1019_norm.vcf,CAD1020:CAD1020_norm.vcf,CAD1021:CAD1021_norm.vcf,CAD1022:CAD1022_norm.vcf,CAD1023:CAD1023_norm.vcf,CAD1025:CAD1025_norm.vcf,CAD1026:CAD1026_norm.vcf,CAD1027:CAD1027_norm.vcf,CAD1028:CAD1028_norm.vcf,CAD1029:CAD1029_norm.vcf,CAD1030:CAD1030_norm.vcf,CAD1031:CAD1031_norm.vcf,CAD1032:CAD1032_norm.vcf,CAD1033:CAD1033_norm.vcf,CAD1034:CAD1034_norm.vcf,CAD1035:CAD1035_norm.vcf,CAD1036:CAD1036_norm.vcf,CAD1037:CAD1037_norm.vcf,CAD1038:CAD1038_norm.vcf,CAD1039:CAD1039_norm.vcf,CAD1040:CAD1040_norm.vcf,CAD1041:CAD1041_norm.vcf,CAD1042:CAD1042_norm.vcf,CAD1043:CAD1043_norm.vcf,CAD1044:CAD1044_norm.vcf,CAD1045:CAD1045_norm.vcf,CAD1046:CAD1046_norm.vcf,CAD1047:CAD1047_norm.vcf,CAD1048:CAD1048_norm.vcf,CAD1049:CAD1049_norm.vcf,CAD1050:CAD1050_norm.vcf,CAD1051:CAD1051_norm.vcf,CAD1052:CAD1052_norm.vcf,CAD1053:CAD1053_norm.vcf,CAD1054:CAD1054_norm.vcf,CAD1055:CAD1055_norm.vcf,CAD1056:CAD1056_norm.vcf,CAD1057:CAD1057_norm.vcf,CAD1058:CAD1058_norm.vcf,CAD1059:CAD1059_norm.vcf,CAD1060:CAD1060_norm.vcf,CAD1061:CAD1061_norm.vcf,CAD1062:CAD1062_norm.vcf,CAD1063:CAD1063_norm.vcf,CAD1064:CAD1064_norm.vcf,CAD1065:CAD1065_norm.vcf,CAD1066:CAD1066_norm.vcf,CAD1067:CAD1067_norm.vcf,CAD1068:CAD1068_norm.vcf,CAD1069:CAD1069_norm.vcf,CAD1070:CAD1070_norm.vcf,CAD1071:CAD1071_norm.vcf,CAD1072:CAD1072_norm.vcf,CAD1073:CAD1073_norm.vcf,CAD1074:CAD1074_norm.vcf,CAD1075:CAD1075_norm.vcf,CAD1076:CAD1076_norm.vcf,CAD1077:CAD1077_norm.vcf,CAD1078:CAD1078_norm.vcf,CAD1079:CAD1079_norm.vcf,CAD1080:CAD1080_norm.vcf,CAD1081:CAD1081_norm.vcf,CAD1082:CAD1082_norm.vcf,CAD1083:CAD1083_norm.vcf,CAD1084:CAD1084_norm.vcf,CAD1085:CAD1085_norm.vcf,CAD1086:CAD1086_norm.vcf,CAD1087:CAD1087_norm.vcf,CAD1088:CAD1088_norm.vcf,CAD1089:CAD1089_norm.vcf,CAD1090:CAD1090_norm.vcf,CAD1091:CAD1091_norm.vcf,CAD1092:CAD1092_norm.vcf,CAD1093:CAD1093_norm.vcf,CAD1094:CAD1094_norm.vcf,CAD1095:CAD1095_norm.vcf,CAD1096:CAD1096_norm.vcf,CAD1097:CAD1097_norm.vcf,CAD1098:CAD1098_norm.vcf,CAD1099:CAD1099_norm.vcf,CAD1100:CAD1100_norm.vcf,CAD1101:CAD1101_norm.vcf,CAD1102:CAD1102_norm.vcf,CAD1103:CAD1103_norm.vcf

# Moving back to the parent directory
cd ..

# --- 4: Splitting the multiallelic records introduced by bayesTyperTools combine
$bcftools norm --multiallelics -any -Ov asmvar_filtering/asmvar_combined.vcf > asmvar_filtering/asmvar_split.vcf

# --- 5: Adding the SVTYPE annotation along with the appropriate header line
# DEPENDENCY : scripts/add_svtype.awk
# DEPENDENCY : illumina_sv_calling/asmvar/svtype_header_line.txt
../../scripts/add_svtype.awk asmvar_filtering/asmvar_split.vcf | \
	awk 'BEGIN {OFS = "\t"} /^#/ {print} !/^#/ {sub("ACO=.*;", "ACO=asmvar;"); print}' | \
	$bcftools annotate --header-lines svtype_header_line.txt -Ov - > asmvar_filtering/asmvar_annotated.vcf

# --- 6: Extracting the SVs >= 50 bp from the annotated VCF file and setting the ID to the name of the caller + line number
# DEPENDENCY : scripts/extract_svs_50.awk
../../scripts/extract_svs_50.awk asmvar_filtering/asmvar_annotated.vcf | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$3 = NR ; print}' | \
	$bcftools annotate --set-id "%INFO/ACO\_%INFO/SVTYPE\_%ID" -Ov - > asmvar_filtering/asmvar_svs.vcf

