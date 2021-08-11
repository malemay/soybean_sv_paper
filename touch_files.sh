# Creating files for sequencing data
mkdir -p illumina_data/raw_fastq

for i in $(cat utilities/all_lines.txt)
do 
	touch illumina_data/raw_fastq/${i}_R1.fastq.gz
	touch illumina_data/raw_fastq/${i}_R2.fastq.gz
done

for i in $(cat utilities/flowcell_names.txt)
do
	touch nanopore_data/${i}.fastq.gz
done

# Creating the reference data files
touch refgenome/bt_decoy_sequences.fasta
touch refgenome/Gmax_508_v4.0.fa
touch refgenome/Gmax_508_v4.0_mit_chlp.fasta
touch refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
touch refgenome/Gmax_508_Wm82.a4.v1.gene_exons.gff3
touch refgenome/Gmax_508_Wm82.a4.v1.gene.gff3
touch refgenome/Gmax_508_Wm82.a4.v1.repeatmasked_assembly_v4.0.gff3

# Create the transposable element data
mkdir -p te_analysis/te_database
touch te_analysis/te_database/SoyBase_TE_Fasta.txt
touch te_analysis/tian2012_tes.txt

# Creating directories that would normally be created by code
mkdir -p illumina_sv_calling/asmvar/asmvar_filtering
mkdir -p sv_genotyping/illumina_svs/sveval_benchmarks/nogeno_RData
mkdir -p sv_genotyping/nanopore_svs/sveval_benchmarks/nogeno_RData
mkdir -p sv_genotyping/illumina_svs/sveval_benchmarks/norepeat_RData
mkdir -p sv_genotyping/illumina_svs/sveval_benchmarks/frequency_RData
mkdir -p sv_genotyping/nanopore_svs/sveval_benchmarks/norepeat_RData
mkdir -p sv_genotyping/combined_svs/sveval_benchmarks/nogeno_RData
mkdir -p sv_genotyping/combined_svs/sveval_benchmarks/norepeat_RData
mkdir -p breakpoint_refinement_analysis/raw_svs/sveval_benchmarks/nogeno_RData/sveval_nogeno_rates.RData
mkdir -p breakpoint_refinement_analysis/refined_svs/sveval_benchmarks/nogeno_RData

