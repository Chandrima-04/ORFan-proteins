SHELL:=/bin/bash

reads_R1:=../../amplicon2/2_qf/4_unmix/d9539_amplicon2_qf_R1.fq.gz
reads_R2:=../../amplicon2/2_qf/4_unmix/d9539_amplicon2_qf_R2.fq.gz
reads_single:=../../amplicon2/2_qf/4_unmix/d9539_amplicon2_qf_single.fq.gz

SAMPLE:= d9539_asm_v1.0_amplicon2
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
reads/$(SAMPLE)_mapped_pe.sam: $(reads_R1) $(reads_R2)| bwa_idx/iva_ctg7_clean1_extented2.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 40 $(basename $|) $^ | samtools sort -O sam -T sort_tmp -@ $(threads) -o $@ -

reads/$(SAMPLE)_mapped_se.sam: $(reads_single) | bwa_idx/iva_ctg7_clean1_extented2.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 40 $(basename $|) $^ | samtools sort -O sam -T sort_tmp -@ $(threads) -o $@ -

$(SAMPLE)_mapped.bam: reads/$(SAMPLE)_mapped_pe.sam reads/$(SAMPLE)_mapped_se.sam
	samtools merge $@ $^

%.bam.bai: %.bam
	samtools index $^
