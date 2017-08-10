SHELL:=/bin/bash

asm_fasta := ../d9539_asm_v1.2.fa

SAMPLE:= d9539_asm_v1.2
threads := 16

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

#all: orf/$(SAMPLE)_emboss_orf0.fa
all: 1_raw_orfs/$(SAMPLE)_emboss_orf1.fa

#*********************************************************
# Predict ORFs with EMBOSS getorf
#*********************************************************
#Gets orfs between two stop codons
# orf/$(SAMPLE)_emboss_orf0.fa: $(asm_fasta)
# 	mkdir -p $(dir $@)
# 	getorf -sequence $< -outseq $@ -circular Y -find 0

#Detects orfs from start to stop codons
1_raw_orfs/$(SAMPLE)_orf.fa: $(asm_fasta)
	mkdir -p $(dir $@)
	getorf -sequence $< -outseq $@ -circular Y -find 1
