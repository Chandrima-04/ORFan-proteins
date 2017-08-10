SHELL:=/bin/bash

SAMPLE:= d9539_asm_v1.2
threads := 16

PfamA_hmm := /scratch/maubar/hmmerdb/Pfam-A/Pfam-A.hmm
vFamA_hmm := /scratch/maubar/hmmerdb/vFam-A/vFam-A_2014.hmm

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: all

all: 2_hmmscan/$(SAMPLE)_orf_hmmscan_PfamA.out
all: 2_hmmscan/$(SAMPLE)_orf_hmmscan_vFamA.out

#*********************************************************
# Run hmmscan Pfam-A against the orfs
#*********************************************************
hmmscan_out_opt= --tblout $(basename $@).tbl --domtblout $(basename $@).domtbl

2_hmmscan/%_hmmscan_PfamA.out: 1_orf/%.fa | $(PfamA_hmm)
	mkdir -p $(dir $@)
	hmmscan --cpu 16 $(hmmscan_out_opt) -o $@ $| $<

2_hmmscan/%_hmmscan_vFamA.out: 1_orf/%.fa | $(vFamA_hmm)
	mkdir -p $(dir $@)
	hmmscan --cpu 16 $(hmmscan_out_opt) -o $@ $| $<
