#!/bin/bash
set -euo pipefail -o verbose

module load ncbi-blast/2.4.0+

OUT_PREFIX=poophage_crispr

blastn -db ../blastdb/crisprdb_170102 -task blastn-short -outfmt 11 -num_threads 8 -query d9539_asm_v1.2.fa -out ${OUT_PREFIX}.asn

blast_formatter -archive ${OUT_PREFIX}.asn -html -out ${OUT_PREFIX}.html

blast_formatter -archive ${OUT_PREFIX}.asn -outfmt '6 qseqid sallseqid pident length mismatch gapopen slen qstart qend sstart send evalue bitscore staxids' -out ${OUT_PREFIX}.tsv
