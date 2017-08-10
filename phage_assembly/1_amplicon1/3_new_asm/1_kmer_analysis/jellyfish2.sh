#!/bin/bash
set -euo pipefail

sample_name=d9539_amplicon1

for kmer_size in 17 31;
do
	jellyfish2 count -s 8G -C -m ${kmer_size} -t 16 -o ${sample_name}_${kmer_size}_raw.count <(zcat ../reads/*.fastq.gz)
	jellyfish2 histo -t 16 ${sample_name}_${kmer_size}_raw.count -o ${sample_name}_k${kmer_size}_raw.hist
	Rscript ~/workspace/meta_illumina_pipeline/scripts/plot_kmer_histogram.R ${sample_name}_k${kmer_size}_raw.hist ${sample_name}_${kmer_size}_raw.pdf
done
