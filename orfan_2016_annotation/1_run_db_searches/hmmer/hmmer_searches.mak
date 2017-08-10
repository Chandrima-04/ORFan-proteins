SHELL := /bin/bash

ifndef threads
threads := 8
endif

.DELETE_ON_ERROR:

.SECONDARY:

.PHONY: nr env_nr

nr: hmmer_vs_nr/$(sample_name)_hmmsearch_vs_nr.txt
env_nr: hmmer_vs_env_nr/$(sample_name)_hmmsearch_vs_env_nr.txt
hmp_pep: hmmer_vs_hmp_pep/$(sample_name)_hmmsearch_vs_hmp_pep.txt
metahit_pep: hmmer_vs_metahit_pep/$(sample_name)_hmmsearch_vs_metahit_pep.txt

#Databases
fasta_nr:= /scratch/maubar/fasta/nr/nr
fasta_env_nr:= /scratch/maubar/fasta/env_nr/env_nr
fasta_hmp_pep:= /scratch/maubar/fasta/hmp/all_pep.fsa
fasta_metahit_pep:= /scratch/maubar/fasta/metahit/UniGene.pep

#***************************************************
# Create profiles from the multiple seq alignments
# for hmmsearch
#***************************************************
hmm_profiles/%.hmm : all_aa_afa/%.fasta
	mkdir -p $(dir $@)
	hmmbuild -n $(sample_name) $@ $<

#***************************************************
# HMMsearch options
#***************************************************
hmmsearch_opts := --cpu $(threads) -E 1 --domE 1 --incE 0.01 --incdomE 0.03 --notextw

hmmer_vs_nr/%_hmmsearch_vs_nr.txt: hmm_profiles/%_AA.hmm
	mkdir -p $(dir $@)
	hmmsearch  $(hmmsearch_opts) -o $@ --tblout $(basename $@).tsv --domtblout $(basename $@).dom $< $(fasta_nr)

hmmer_vs_env_nr/%_hmmsearch_vs_env_nr.txt: hmm_profiles/%_AA.hmm
	mkdir -p $(dir $@)
	hmmsearch  $(hmmsearch_opts) -o $@ --tblout $(basename $@).tsv --domtblout $(basename $@).dom $< $(fasta_env_nr)

hmmer_vs_hmp_pep/%_hmmsearch_vs_hmp_pep.txt: hmm_profiles/%_AA.hmm
	mkdir -p $(dir $@)
	hmmsearch  $(hmmsearch_opts) -o $@ --tblout $(basename $@).tsv --domtblout $(basename $@).dom $< $(fasta_hmp_pep)

hmmer_vs_metahit_pep/%_hmmsearch_vs_metahit_pep.txt: hmm_profiles/%_AA.hmm
	mkdir -p $(dir $@)
	hmmsearch  $(hmmsearch_opts) -o $@ --tblout $(basename $@).tsv --domtblout $(basename $@).dom $< $(fasta_metahit_pep)
