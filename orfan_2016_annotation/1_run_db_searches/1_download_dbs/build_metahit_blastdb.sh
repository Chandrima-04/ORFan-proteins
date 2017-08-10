#!/bin/bash
mkdir -p blastdb/metahit/{cds,pep}

makeblastdb -input_type fasta -out blastdb/metahit/cds/metahit_cds -title metahit_cds -dbtype nucl -in fasta/metahit/UniGene.cds
makeblastdb -input_type fasta -out blastdb/metahit/pep/metahit_pep -title metahit_pep -dbtype prot -in fasta/metahit/UniGene.pep
