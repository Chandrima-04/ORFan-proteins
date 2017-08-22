#!/bin/bash
set -euo pipefail -o verbose
IFS=$'\n\t'

module load ncbi-blast/2.4.0+
module load samtools/1.3

#1. Create blast index of viruses Indexes
virus_dir=0_virus_fasta
mkdir -p ${virus_dir}/blastdb
cat ${virus_dir}/*.fasta > ${virus_dir}/blastdb/viruses.fasta

makeblastdb -input_type fasta -out ${virus_dir}/blastdb/viruses -title viruses -dbtype nucl -in ${virus_dir}/blastdb/viruses.fasta

blastdb=${virus_dir}/blastdb/viruses

#Run blastn
mkdir -p 1_blastn
cd 1_blastn

in_fasta=../../0_data/orfans/ALL_CLUSTERS.fa
out_prefix=orfan_to_viruses

make -r -f ../blastn.mak in_fasta=../${in_fasta} blastdb=../${blastdb} out_prefix=${out_prefix} tsv

#Query coverage should be at least 30%
awk '{if($5 > 30) print $0}' tsv/${out_prefix}.tsv > tsv/${out_prefix}_filt.tsv

mkdir -p plots/
python ../plot_blastn_hits.py tsv/${out_prefix}_filt.tsv tsv/blast_tsv_columns.txt plots/
