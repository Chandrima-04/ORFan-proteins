#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

module load hmmer/3.1b2

IN_FASTA=454_seqs_phage.fa
HMM_DB=/labcommon/db/hmmerdb/Pfam-A_r31/Pfam-A.hmm
OUT_PREFIX=454_seqs_phage_PfamA

hmmscan --cpu 16 --noali --tblout ${OUT_PREFIX}.tbl \
	--domtblout ${OUT_PREFIX}.domtbl \
	-o ${OUT_PREFIX}.out ${HMM_DB} ${IN_FASTA}
