#!/bin/bash
for x in $(ls all_nt_unaligned);
do
echo "Processing ${x%_DNA.fa}"
make -r -f blast_searches.mak sample_name=${x%_DNA.fa} env_nt env_nr
make -r -f blast_searches.mak sample_name=${x%_DNA.fa} nr
make -r -f blast_searches.mak sample_name=${x%_DNA.fa} nt
make -r -f blast_searches.mak sample_name=${x%_DNA.fa} hmp_nuc
make -r -f blast_searches.mak sample_name=${x%_DNA.fa} hmp_pep
make -r -f blast_searches.mak sample_name=${x%_DNA.fa} metahit_cds metahit_pep
make -r -f hmmer_searches.mak sample_name=${x%_DNA.fa} env_nr
make -r -f hmmer_searches.mak sample_name=${x%_DNA.fa} nr
make -r -f hmmer_searches.mak sample_name=${x%_DNA.fa} hmp_pep
make -r -f hmmer_searches.mak sample_name=${x%_DNA.fa} metahit_pep
echo
done
