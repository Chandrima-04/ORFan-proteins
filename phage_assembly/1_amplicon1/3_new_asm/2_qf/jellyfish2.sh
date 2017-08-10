#!/bin/bash
set -euo pipefail

sample_name=$1
read_folder=$2

for kmer_size in 17 31;
do
	jellyfish2 count -s 8G -C -m ${kmer_size} -t 16 -o ${sample_name}_k${kmer_size}.count <(zcat ${read_folder}/*.f*q.gz)
	jellyfish2 histo -t 16 ${sample_name}_k${kmer_size}.count -o ${sample_name}_k${kmer_size}.hist
	Rscript ~/workspace/meta_illumina_pipeline/scripts/plot_kmer_histogram.R ${sample_name}_k${kmer_size}.hist ${sample_name}_k${kmer_size}.pdf
done
