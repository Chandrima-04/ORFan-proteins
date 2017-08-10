SHELL:=/bin/bash

reads_pe:=../../amplicon1/3_new_asm/2_qf/4_phix_rm/d9539_amplicon1_qf_pe.fq.gz
reads_se:=../../amplicon1/3_new_asm/2_qf/4_phix_rm/d9539_amplicon1_qf_se.fq.gz

SAMPLE:= d9539_asm_v1.0_amplicon1
threads := 16

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: ../contig/iva_ctg7_clean1_extended2.fa.fai $(SAMPLE)_mapped.bam.bai

#*********************************************************
# Create a fasta index from the contigs
#*********************************************************
%.fa.fai: %.fa
	samtools faidx $^

#*************************************************************
# Map reads back to contig
#*************************************************************
#Create bwa idx from contig
bwa_idx/iva_ctg7_clean1_extented2.bwt : ../contig/iva_ctg7_clean1_extended2.fa
	mkdir -p $(dir $@)
	bwa index -p $(basename $@) $^

#Map reads against contigs
reads/$(SAMPLE)_mapped_pe.sam: $(reads_pe) | bwa_idx/iva_ctg7_clean1_extented2.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 40 -p $(basename $|) $< | samtools sort -O sam -T sort_tmp -@ $(threads) -o $@ -

reads/$(SAMPLE)_mapped_se.sam: $(reads_se) | bwa_idx/iva_ctg7_clean1_extented2.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 40 $(basename $|) $< | samtools sort -O sam -T sort_tmp -@ $(threads) -o $@ -

$(SAMPLE)_mapped.bam: reads/$(SAMPLE)_mapped_pe.sam reads/$(SAMPLE)_mapped_se.sam
	samtools merge $@ $^

%.bam.bai: %.bam
	samtools index $^
