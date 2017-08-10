SHELL := /bin/bash

ifndef threads
threads := 8
endif

.DELETE_ON_ERROR:

.PHONY: nt nr env_nt env_nr hmp_nuc hmp_pep metahit_cds metahit_pep

nt: $(addprefix 2_blast/$(sample_name)_blastn_vs_nt,.html .tsv)
env_nt: $(addprefix 2_blast/$(sample_name)_blastn_vs_env_nt,.html .tsv)
hmp_nuc: $(addprefix 2_blast/$(sample_name)_blastn_vs_hmp_nuc,.html .tsv)
metahit_cds: $(addprefix 2_blast/$(sample_name)_blastn_vs_metahit_cds,.html .tsv)

nr: $(addprefix 2_blast/$(sample_name)_blastp_vs_nr,.html .tsv)
env_nr: $(addprefix 2_blast/$(sample_name)_blastp_vs_env_nr,.html .tsv)
hmp_pep: $(addprefix 2_blast/$(sample_name)_blastp_vs_hmp_pep,.html .tsv)
metahit_pep: $(addprefix 2_blast/$(sample_name)_blastp_vs_metahit_pep,.html .tsv)

#Databases
blastdb_nt:= /scratch/maubar/blastdb/nt/nt
blastdb_env_nt:= /scratch/maubar/blastdb/env_nt/env_nt
blastdb_hmp_nuc:= /scratch/maubar/blastdb/hmp/nuc/hmp_nuc
blastdb_metahit_cds:= /scratch/maubar/blastdb/metahit/cds/metahit_cds

blastdb_nr:= /scratch/maubar/blastdb/nr/nr
blastdb_env_nr:= /scratch/maubar/blastdb/env_nr/env_nr
blastdb_hmp_pep:= /scratch/maubar/blastdb/hmp/pep/hmp_pep
blastdb_metahit_pep:= /scratch/maubar/blastdb/metahit/pep/metahit_pep

#Blast parameters
blast_params:= -num_threads $(threads) -max_target_seqs 10 -outfmt 11

#***************************************************
# BLASTn searches
#***************************************************
2_blast/%_blastn_vs_nt.asn: 1_orf/%.fa
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb_nt) -query $< -out $@

2_blast/%_blastn_vs_env_nt.asn: 1_orf/%.fa
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb_env_nt) -query $< -out $@

2_blast/%_blastn_vs_hmp_nuc.asn: 1_orf/%.fa
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb_hmp_nuc) -query $< -out $@

2_blast/%_blastn_vs_metahit_cds.asn: 1_orf/%.fa
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb_metahit_cds) -query $< -out $@

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

#********
# Formatters
#**********
%.html: %.asn
	mkdir -p $(dir $@)
	blast_formatter -archive $< -html -out $@

%.tsv: %.asn
	mkdir -p $(dir $@)
	blast_formatter -archive $< -outfmt 6 -out $@
