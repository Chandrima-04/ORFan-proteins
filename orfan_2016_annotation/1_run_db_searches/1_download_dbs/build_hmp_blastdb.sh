#!/bin/bash
mkdir -p blastdb/hmp/{nuc,pep}

makeblastdb -input_type fasta -out blastdb/hmp/nuc/hmp_nuc -title hmp_nuc -dbtype nucl -in fasta/hmp/all_nuc.fsa
makeblastdb -input_type fasta -out blastdb/hmp/pep/hmp_pep -title hmp_pep -dbtype prot -in fasta/hmp/all_pep.fsa
