SHELL:=/bin/bash

amplicon1_reads:=$(wildcard reads/amplicon1/*_1.fastq.gz) $(wildcard reads/amplicon1/*_2.fastq.gz)
amplicon2_reads:=$(wildcard reads/amplicon2/*_R1_*.fastq.gz) $(wildcard reads/amplicon2/*_R2_*.fastq.gz)

SAMPLE:= d9539_asm_v1.2
threads := 16

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

.PHONY: indexes all

all: indexes bamfiles

indexes: idx/d9539_asm_v1.2.prj idx/d9539_asm_v1.2.bwt
bam_files: mapped/$(SAMPLE)_amplicon1.bam.bai
bam_files: mapped/$(SAMPLE)_amplicon2.bam.bai
bam_files: mapped/$(SAMPLE)_cluster179.bam.bai
bam_files: mapped/$(SAMPLE)_cluster179a_DNA.bam.bai
bam_files: mapped/$(SAMPLE)_cluster179b_DNA.bam.bai
bam_files: mapped/$(SAMPLE)_env_nt.bam.bai
bam_files: mapped/$(SAMPLE)_metahit.bam.bai


#*********************************************************
# Create a fasta index from the contigs
#*********************************************************
%.fa.fai: %.fa
	samtools faidx $^

#*********************************************************
# Create mapper indexes
#*********************************************************
idx/%.prj: %.fa
	mkdir -p $(dir $@)
	lastdb -cR01 $(basename $@) $<

idx/%.bwt: %.fa
	mkdir -p $(dir $@)
	bwa index -p $(basename $@) $<


#*************************************************************
# Map reads back to contig
#*************************************************************
#Map reads against contigs
mapped/$(SAMPLE)_amplicon1.bam: $(amplicon1_reads) | idx/d9539_asm_v1.2.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 40 $(basename $|) $^ | samtools view -Sb - | samtools sort -O bam -T sort_tmp -@ 16 -o $@ -

mapped/$(SAMPLE)_amplicon2.bam: $(amplicon2_reads) | idx/d9539_asm_v1.2.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 40 $(basename $|) $^ | samtools view -Sb - | samtools sort -O bam -T sort_tmp -@ 16 -o $@ -


mapped/$(SAMPLE)_%.maf: reads/cluster179/%.fa | idx/d9539_asm_v1.2.prj
	mkdir -p $(dir $@)
	parallel-fasta -j 8 "lastal $(basename $| )" < $< > $@

mapped/$(SAMPLE)_%.maf: reads/homologs/d9539_%_homologs.fasta | idx/d9539_asm_v1.2.prj
	mkdir -p $(dir $@)
	parallel-fasta -j 8 "lastal $(basename $| )" < $< > $@

%.bam: %.maf
	maf-convert sam $< | cat sam.header - | samtools view -hSb - | samtools sort -O bam -T sort_tmp -@ $(threads) -o $@ -

%.bam.bai: %.bam
	samtools index $^

