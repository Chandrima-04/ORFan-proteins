#!/bin/bash
set -euo verbose

for x in $(ls all_nt_unaligned);
do
echo "Processing ${x%_DNA.fa}"
make -r -f blast_searches.mak sample_name=${x%_DNA.fa} blastn_igc_cds
make -r -f blast_searches.mak sample_name=${x%_DNA.fa} blastp_igc_pep
make -r -f blast_searches.mak sample_name=${x%_DNA.fa} hmmer_igc_pep
echo
done
