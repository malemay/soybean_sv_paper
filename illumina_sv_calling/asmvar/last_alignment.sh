
# Creating a few variables to make scripting easier
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
ref=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta
lastal=$1
lastsplit=$2

# Looping over the samples
lines=$(cat ../../utilities/all_lines.txt)

# DEPENDENCY : SOAPdenovo 2 assemblies
for i in $lines
do
    # I use option -m 0.01 to keep the previous default setting of last-926
    $lastal -D1000 -P2 -Q0 -e20 -j4 -v $ref soap_assembly/${i}/${i}.contig | $lastsplit -m 0.01 -s30 -v > last_alignment/${i}_soapdenovo2_fm.maf
done

