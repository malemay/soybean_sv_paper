#!/bin/bash

# Getting the path to the manta executable from the command line
manta=$1

# The analysis will be run in 10 batches of 10/11 samples
# Which samples are part of what batch has been determined randomly

# DEPENDENCY : Illumina reads aligned with BWA
# DEPENDENCY : illumina_sv_calling/manta/manta_config.txt

# Processing the first batch
mkdir -p manta_sample1
cd manta_sample1

# Running the actual command; output will be written to MantaWorkflow as per the default settings
$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1033/CAD1033_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1085/CAD1085_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1021/CAD1021_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1077/CAD1077_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1030/CAD1030_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1067/CAD1067_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1045/CAD1045_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1017/CAD1017_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1012/CAD1012_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1044/CAD1044_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j 2 -g 4
cd ../..

# Processing the second batch
mkdir -p manta_sample2
cd manta_sample2

# Running the actual command; output will be written to MantaWorkflow as per the default settings
$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1099/CAD1099_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1048/CAD1048_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1097/CAD1097_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1094/CAD1094_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1065/CAD1065_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1035/CAD1035_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1013/CAD1013_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1082/CAD1082_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1081/CAD1081_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1001/CAD1001_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j 2 -g 4
cd ../..

# Processing the third batch
mkdir -p manta_sample3
cd manta_sample3

# Running the actual command; output will be written to MantaWorkflow as per the default settings
$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1091/CAD1091_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1086/CAD1086_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1020/CAD1020_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1075/CAD1075_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1038/CAD1038_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1023/CAD1023_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1002/CAD1002_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1026/CAD1026_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1029/CAD1029_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1096/CAD1096_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j 4 -g 4
cd ../..

# Processing the fourth batch
mkdir -p manta_sample4
cd manta_sample4

# Running the actual command; output will be written to MantaWorkflow as per the default settings
$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1100/CAD1100_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1059/CAD1059_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1061/CAD1061_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1019/CAD1019_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1014/CAD1014_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1007/CAD1007_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1052/CAD1052_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1055/CAD1055_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1071/CAD1071_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1040/CAD1040_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j 2 -g 4
cd ../..

# Processing the fifth batch
mkdir -p manta_sample5
cd manta_sample5

# Running the actual command; output will be written to MantaWorkflow as per the default settings
$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1056/CAD1056_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1054/CAD1054_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1072/CAD1072_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1039/CAD1039_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1060/CAD1060_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1006/CAD1006_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1047/CAD1047_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1036/CAD1036_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1025/CAD1025_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1093/CAD1093_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j 4 -g 4
cd ../..

# Processing the sixth batch
mkdir -p manta_sample6
cd manta_sample6

# Running the actual command; output will be written to MantaWorkflow as per the default settings
$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1087/CAD1087_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1049/CAD1049_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1008/CAD1008_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1068/CAD1068_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1074/CAD1074_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1078/CAD1078_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1073/CAD1073_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1043/CAD1043_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1042/CAD1042_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1046/CAD1046_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j 2 -g 4
cd ../..

# Processing the seventh batch
mkdir -p manta_sample7
cd manta_sample7

$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1010/CAD1010_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1022/CAD1022_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1090/CAD1090_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1084/CAD1084_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1069/CAD1069_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1092/CAD1092_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1080/CAD1080_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1098/CAD1098_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1095/CAD1095_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1064/CAD1064_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j 4 -g 4
cd ../..

# Processing the eighth batch
mkdir -p manta_sample8
cd manta_sample8

$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1089/CAD1089_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1051/CAD1051_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1070/CAD1070_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1031/CAD1031_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1018/CAD1018_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1027/CAD1027_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1101/CAD1101_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1058/CAD1058_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1028/CAD1028_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1062/CAD1062_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j 4 -g 4
cd ../..

# Processing the ninth batch
mkdir -p manta_sample9
cd manta_sample9

$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1079/CAD1079_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1011/CAD1011_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1034/CAD1034_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1005/CAD1005_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1032/CAD1032_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1057/CAD1057_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1009/CAD1009_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1102/CAD1102_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1103/CAD1103_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1015/CAD1015_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1004/CAD1004_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j 4 -g 4
cd ../..

# Processing the tenth batch
mkdir -p manta_sample10
cd manta_sample10

$manta --config ../manta_config.txt --bam ../../../illumina_data/aligned_reads/CAD1066/CAD1066_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1083/CAD1083_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1037/CAD1037_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1016/CAD1016_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1041/CAD1041_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1050/CAD1050_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1053/CAD1053_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1063/CAD1063_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1003/CAD1003_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1088/CAD1088_all.sort.bam \
    --bam ../../../illumina_data/aligned_reads/CAD1076/CAD1076_all.sort.bam

# Now launching the analysis using the executable that has been created
cd MantaWorkflow
./runWorkflow.py -j8
cd ../..

