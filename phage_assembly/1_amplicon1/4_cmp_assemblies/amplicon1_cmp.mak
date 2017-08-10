SHELL:=/bin/bash

sample_name:= d9539_amplicon2

paired_end:=../3_new_asm/4_phix_rm/$(sample_name)_qf_pe.fq.gz
single:=../3_new_asm/4_phix_rm/$(sample_name)_qf_se.fq.gz

threads := 16

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: $(addsuffix .bam.bai,ctgs/megahit ctgs/fermi )

#*********************************************************
# Get contig.00007 from the IVA assembly which is
#    the one we are interested in
#*********************************************************
asm1_iva_ctg7/iva_ctg7.fa: ../1_asm/iva-raw/contigs.fasta
	mkdir -p $(dir $@)
	seqtk subseq $< <(echo "contig.00007") > $@

#*********************************************************
# Create a last indexes from the Contig 7
#*********************************************************
asm1_iva_ctg7/last_idx/iva_ctg7.prj: asm1_iva_ctg7/iva_ctg7.fa
	mkdir -p $(dir $@)
	lastdb -cR01 $(basename $@) $^

#*************************************************************
# Align amplicon2 contigs against amplicon1 contig 7
#*************************************************************
ctgs/megahit.maf: ../3_new_asm/3_asm/megahit/final.contigs.fa | asm1_iva_ctg7/last_idx/iva_ctg7.prj
ctgs/fermi.maf:   ../3_new_asm/3_asm/fermi/fermi_contigs.fa   | asm1_iva_ctg7/last_idx/iva_ctg7.prj
# ctgs/spades.maf:  ../3_new_asm/spades/contigs.fasta     | amplicon1_iva_ctg7/last_idx/iva_ctg7.prj

%.maf:
	mkdir -p $(dir $@)
	parallel-fasta -j 8 "lastal $(basename $| )" < <(seqtk seq -L 500 $< ) > $@

%.bam: %.maf
	maf-convert sam $< | cat sam.header - | samtools view -hSb - | samtools sort -O bam -T sort_tmp -@ $(threads) -o $@ -

%.bam.bai : %.bam
	samtools index $^

#*************************************************************
# Map reads back to assembly
#*************************************************************
#Map reads against contigs
reads/mapped_pe.sam: $(paired_end) | amplicon1_iva_ctg7/bwa_idx/iva_ctg7.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 30 $(basename $|) <(seqtk seq -1 $<) <(seqtk seq -2 $<) | samtools sort -O sam -T sort_tmp -@ $(threads) -o $@ -

reads/mapped_se.sam: $(single) | amplicon1_iva_ctg7/bwa_idx/iva_ctg7.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 30 $(basename $|) $^ | samtools sort -O sam -T sort_tmp -@ $(threads) -o $@ -

reads/mapped.bam: reads/mapped_pe.sam reads/mapped_se.sam
	samtools merge $@ $^

%.bam.bai: %.bam
	samtools index $^
