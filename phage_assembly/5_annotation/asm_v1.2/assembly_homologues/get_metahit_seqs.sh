#!/bin/bash

alias extract_first_fasta_seq='awk '"\'"'BEGIN{seq=0} $0~/^>/{if(seq == 0){seq=1; print $0} else exit} $0 !~/^>/{print $0}'"\'"

parallel grep -F -A 10 {} fasta/metahit/UniGene.cds :::: d9539_metahit_homologs.txt


xargs --arg-file d9539_metahit_homologs.txt -I{} grep -F -A 1 {} fasta/metahit/UniGene.cds > d9539_metahit_homologs.fasta