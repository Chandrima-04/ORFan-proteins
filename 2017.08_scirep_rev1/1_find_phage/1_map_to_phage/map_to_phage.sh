#!/bin/bash
set -euo pipefail

module load bwa/0.7.15
module load samtools/1.3

DATA_DIR=../../0_data/

BWA_IDX=${DATA_DIR}/bwa_idx/d9539_asm_v1.2 
IN_FASTA=${DATA_DIR}/454_seqs/myunknown_goodclean5_prinsgood.fa

bwa mem -T30 -M -t 16 ${BWA_IDX} ${IN_FASTA} | samtools view -hSb -F4 - | samtools sort - > 454_seqs_to_phage.bam

samtools index *.bam
