SHELL:=/bin/bash

reads_pe:=../3_new_asm/2_qf/4_phix_rm/*_pe.fq.gz
reads_se:=../3_new_asm/2_qf/4_phix_rm/*_se.fq.gz

SAMPLE:= d9539_amplicon1_scf53612
threads := 16

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: ctgs/masurca_scf53612.fa.fai reads/mapped.bam.bai

#*********************************************************
# Get contig.00007 from the IVA assembly which is
#    the one we are interested in
#*********************************************************
ctgs/masurca_scf53612.fa: ../1_asm/masurca-carlos/A7VRF_P1141.masurca-K71.4Kb.fasta
	mkdir -p $(dir $@)
	seqtk subseq $< <(echo "scf7180000053612") > $@

#*********************************************************
# Create a fasta index from the contigs
#*********************************************************
%.fa.fai: %.fa
	samtools faidx $^

#*************************************************************
# Map reads back to contig
#*************************************************************
#Create bwa idx from contig
bwa_idx/masurca_scf53612.bwt : ctgs/masurca_scf53612.fa
	mkdir -p $(dir $@)
	bwa index -p $(basename $@) $^

#Map reads against contigs
reads/mapped_pe.sam: $(reads_pe) | bwa_idx/masurca_scf53612.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 40 $(basename $|) <(seqtk seq -1 $<) <(seqtk seq -2 $<) | samtools sort -O sam -T sort_tmp -@ $(threads) -o $@ -

reads/mapped_se.sam: $(reads_se) | bwa_idx/masurca_scf53612.bwt
	mkdir -p $(dir $@)
	bwa mem -t $(threads) -M -T 40 $(basename $|) $^ | samtools sort -O sam -T sort_tmp -@ $(threads) -o $@ -

reads/mapped.bam: reads/mapped_pe.sam reads/mapped_se.sam
	samtools merge $@ $^

%.bam.bai: %.bam
	samtools index $^
