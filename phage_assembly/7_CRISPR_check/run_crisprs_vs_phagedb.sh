#!/bin/bash
set -euo pipefail -o verbose

module load ncbi-blast/2.4.0+

OUT_PREFIX=crisprs_vs_poophage

blastn -db crisprs_vs_phagedb/poophage_v1.2 -task blastn-short -soft_masking false -evalue 0.5 -outfmt 11 -num_threads 8 -query ../format_crisprdb/CRISPRdb_spacer_170102.fa -out ${OUT_PREFIX}.asn

blast_formatter -archive ${OUT_PREFIX}.asn -html -out ${OUT_PREFIX}.html

blast_formatter -archive ${OUT_PREFIX}.asn -outfmt '6 qseqid sallseqid pident length mismatch gapopen qlen qcovs slen qstart qend sstart send evalue bitscore' -out ${OUT_PREFIX}.tsv
