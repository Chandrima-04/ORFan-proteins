#!/bin/bash
set -euo pipefail

module load bowtie2/2.3.0
module load metaphlan2/2.6.0

mpa2_pkl=${mpa2_dir}/db_v20/mpa_v20_m200.pkl
mpa2_bowtie2db=${mpa2_dir}/db_v20/mpa_v20_m200

IN_FASTA=../../0_data/454_seqs/myunknown_goodclean5_prinsgood.fa
OUT_PREFIX=454_reads_mpa2

metaphlan2.py --mpa_pkl ${mpa2_pkl} --bowtie2db ${mpa2_bowtie2db} \
	--nproc 16 --input_type multifasta --sample_id ${OUT_PREFIX} \
	--bowtie2out ${OUT_PREFIX}.bowtie2.bz2 -s ${OUT_PREFIX}.sam \
	 ${IN_FASTA} > ${OUT_PREFIX}.txt
	
