SHELL := /bin/bash

ifndef threads
threads := 8
endif

.DELETE_ON_ERROR:

.SECONDARY:

.PHONY: nr env_nr hmp_pep metahit_pep

nr: $(addprefix 2_blast/$(sample_name)_blastp_vs_nr,.html .tsv)
env_nr: $(addprefix 2_blast/$(sample_name)_blastp_vs_env_nr,.html .tsv)
hmp_pep: $(addprefix 2_blast/$(sample_name)_blastp_vs_hmp_pep,.html .tsv)
metahit_pep: $(addprefix 2_blast/$(sample_name)_blastp_vs_metahit_pep,.html .tsv)
metahit_2014_pep: $(addprefix 2_blast/$(sample_name)_blastp_vs_metahit_2014_pep,.html .tsv)

blastdb_nr:= /scratch/maubar/blastdb/nr/nr
blastdb_env_nr:= /scratch/maubar/blastdb/env_nr/env_nr
blastdb_hmp_pep:= /scratch/maubar/blastdb/hmp/pep/hmp_pep
blastdb_metahit_pep:= /scratch/maubar/blastdb/metahit/pep/metahit_pep
blastdb_metahit_2014_pep:= /scratch/maubar/blastdb/metahit_2014/igc_pep

#Blast parameters
blast_params:= -num_threads $(threads) -max_target_seqs 10 -outfmt 11

#***************************************************
# BLASTp searches
#***************************************************
2_blast/%_blastp_vs_nr.asn: 1_orf/%.fa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_nr) -query $< -out $@

2_blast/%_blastp_vs_env_nr.asn: 1_orf/%.fa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_env_nr) -query $< -out $@

2_blast/%_blastp_vs_hmp_pep.asn: 1_orf/%.fa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_hmp_pep) -query $< -out $@

2_blast/%_blastp_vs_metahit_pep.asn: 1_orf/%.fa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_metahit_pep) -query $< -out $@

2_blast/%_blastp_vs_metahit_2014_pep.asn: 1_orf/%.fa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_metahit_2014_pep) -query $< -out $@

#******************
# Formatters
#**********
%.html: %.asn
	mkdir -p $(dir $@)
	blast_formatter -archive $< -html -out $@

%.tsv: %.asn
	mkdir -p $(dir $@)
	blast_formatter -archive $< -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sgi staxids sscinames scomnames qcovs stitle' -out $@
