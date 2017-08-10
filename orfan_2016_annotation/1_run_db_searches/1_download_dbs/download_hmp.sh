#!/bin/bash
#Download the references
wget -N http://downloads.hmpdacc.org/data/reference_genomes/all_nuc_20141006.tar.gz
wget -N http://downloads.hmpdacc.org/data/reference_genomes/all_pep_20141006.tar.gz

#Extract the sequences
tar -xzf all_nuc_20141006.tar.gz
tar -xzf all_pep_20141006.tar.gz

#For some reason, the fasta files are nested in a series of folders
# so extract the deepest nested folder with the folder and remove the crap
mv local/db/repository/ncbi/dacc_reference_genomes/20141006/all_nuc_20141006 ./ && rm -r local/
mv local/db/repository/ncbi/dacc_reference_genomes/20141006/all_pep_20141006 ./ && rm -r local/

#Concatenate everything into a single fastas to create the blastdb
cat all_nuc_20141006/*.fsa > all_nuc.fsa
cat all_pep_20141006/*.fsa > all_pep.fsa
