SHELL := /bin/bash

ifndef threads
threads := 8
endif

.DELETE_ON_ERROR:

.SECONDARY:

.PHONY: nt nr env_nt env_nr hmp_nuc hmp_pep metahit_cds metahit_pep

nt: $(addprefix blastn_vs_nt/$(sample_name)_blastn_vs_nt,.html .tsv)
env_nt: $(addprefix blastn_vs_env_nt/$(sample_name)_blastn_vs_env_nt,.html .tsv)
hmp_nuc: $(addprefix blastn_vs_hmp_nuc/$(sample_name)_blastn_vs_hmp_nuc,.html .tsv)
metahit_cds: $(addprefix blastn_vs_metahit_cds/$(sample_name)_blastn_vs_metahit_cds,.html .tsv)

nr: $(addprefix blastp_vs_nr/$(sample_name)_blastp_vs_nr,.html .tsv)
env_nr: $(addprefix blastp_vs_env_nr/$(sample_name)_blastp_vs_env_nr,.html .tsv)
hmp_pep: $(addprefix blastp_vs_hmp_pep/$(sample_name)_blastp_vs_hmp_pep,.html .tsv)
metahit_pep: $(addprefix blastp_vs_metahit_pep/$(sample_name)_blastp_vs_metahit_pep,.html .tsv)

#Databases
blastdb_nt:= /scratch/maubar/blastdb/nt/nt
blastdb_nr:= /scratch/maubar/blastdb/nr/nr

blastdb_env_nt:= /scratch/maubar/blastdb/env_nt/env_nt
blastdb_env_nr:= /scratch/maubar/blastdb/env_nr/env_nr

blastdb_hmp_nuc:= /scratch/maubar/blastdb/hmp/nuc/hmp_nuc
blastdb_hmp_pep:= /scratch/maubar/blastdb/hmp/pep/hmp_pep

blastdb_metahit_cds:= /scratch/maubar/blastdb/metahit/cds/metahit_cds
blastdb_metahit_pep:= /scratch/maubar/blastdb/metahit/pep/metahit_pep

#Blast parameters
blast_params:= -num_threads $(threads) -max_target_seqs 10 -outfmt 11

#***************************************************
# BLASTn searches
#***************************************************
blastn_vs_nt/%_blastn_vs_nt.asn: all_nt_unaligned/%_DNA.fa
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb_nt) -query $< -out $@

blastn_vs_env_nt/%_blastn_vs_env_nt.asn: all_nt_unaligned/%_DNA.fa
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb_env_nt) -query $< -out $@

blastn_vs_hmp_nuc/%_blastn_vs_hmp_nuc.asn: all_nt_unaligned/%_DNA.fa
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb_hmp_nuc) -query $< -out $@

blastn_vs_metahit_cds/%_blastn_vs_metahit_cds.asn: all_nt_unaligned/%_DNA.fa
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb_metahit_cds) -query $< -out $@

#***************************************************
# BLASTp searches
#***************************************************
blastp_vs_nr/%_blastp_vs_nr.asn: all_aa_unaligned/%_AA.fa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_nr) -query $< -out $@

blastp_vs_env_nr/%_blastp_vs_env_nr.asn: all_aa_unaligned/%_AA.fa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_env_nr) -query $< -out $@

blastp_vs_hmp_pep/%_blastp_vs_hmp_pep.asn: all_aa_unaligned/%_AA.fa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_hmp_pep) -query $< -out $@

blastp_vs_metahit_pep/%_blastp_vs_metahit_pep.asn: all_aa_unaligned/%_AA.fa
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
	blast_formatter -archive $< -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sgi staxids sscinames scomnames qcovs stitle' -out $@
