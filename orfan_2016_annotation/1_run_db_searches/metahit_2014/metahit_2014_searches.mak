SHELL := /bin/bash

ifndef threads
threads := 8
endif

.DELETE_ON_ERROR:

.SECONDARY:

.PHONY: blastn_igc_cds blastp_igc_pep hmmer_igc_pep

blastn_igc_cds: $(addprefix blastn_vs_metahit_2014_cds/$(sample_name)_blastn_vs_metahit_2014_cds,.html .tsv)

blastp_igc_pep: $(addprefix blastp_vs_metahit_2014_pep/$(sample_name)_blastp_vs_metahit_2014_pep,.html .tsv)

hmmer_igc_pep: hmmer_vs_metahit_2014_pep/$(sample_name)_hmmsearch_vs_metahit_2014_pep.txt

#Databases

blastdb_metahit_cds:= /scratch/maubar/blastdb/metahit_2014/igc_cds
blastdb_metahit_pep:= /scratch/maubar/blastdb/metahit_2014/igc_pep
fasta_metahit_pep := /scratch/maubar/fasta/metahit_2014/IGC.pep

#Blast parameters
blast_params:= -num_threads $(threads) -max_target_seqs 10 -outfmt 11

#***************************************************
# BLASTn
#***************************************************
blastn_vs_metahit_2014_cds/%_blastn_vs_metahit_2014_cds.asn: all_nt_unaligned/%_DNA.fa
	mkdir -p $(dir $@)
	blastn $(blast_params) -db $(blastdb_metahit_cds) -query $< -out $@

#***************************************************
# BLASTp
#***************************************************
blastp_vs_metahit_2014_pep/%_blastp_vs_metahit_2014_pep.asn: all_aa_unaligned/%_AA.fa
	mkdir -p $(dir $@)
	blastp $(blast_params) -db $(blastdb_metahit_pep) -query $< -out $@

#***************************************************
# HMMsearch
#***************************************************
hmmsearch_opts := --cpu $(threads) -E 1 --domE 1 --incE 0.01 --incdomE 0.03 --notextw

hmmer_vs_metahit_2014_pep/%_hmmsearch_vs_metahit_2014_pep.txt: hmm_profiles/%_AA.hmm
	mkdir -p $(dir $@)
	hmmsearch $(hmmsearch_opts) -o $@ --tblout $(basename $@).tsv --domtblout $(basename $@).dom $< $(fasta_metahit_pep) 


#*********************
# Formatters
#*******************
%.html: %.asn
	mkdir -p $(dir $@)
	blast_formatter -archive $< -html -out $@

%.tsv: %.asn
	mkdir -p $(dir $@)
	blast_formatter -archive $< -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sgi staxids sscinames scomnames qcovs stitle' -out $@
