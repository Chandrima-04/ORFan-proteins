SHELL:=/bin/bash

R1:=$(wildcard ../1_asm/nesoni/*_R1.fq.gz)
R2:=$(wildcard ../1_asm/nesoni/*_R2.fq.gz)
singles:=$(wildcard ../1_asm/nesoni/*_single.fq.gz)

SAMPLE:= d9539_amplicon1
threads := 16

FREEBAYES_BIN := freebayes

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: ctgs/iva_ctg7.fa.fai mapped_sort_rmdup.bam.bai mapped_sort_rmdup.vcf

#*********************************************************
# Get contig.00007 from the IVA assembly which is
#    the one we are interested in
#*********************************************************
ctgs/iva_ctg7.fa: ../1_asm/iva-raw/contigs.fasta
	mkdir -p $(dir $@)
	seqtk subseq $< <(echo "contig.00007") > $@

#*********************************************************
# Create a fasta index from the contigs
#*********************************************************
%.fa.fai: %.fa
	samtools faidx $^

#*************************************************************
# Map reads back to assembly
#*************************************************************
#Create bwa idx from assembly
bwa_idx/iva_ctg7.bwt : ctgs/iva_ctg7.fa
	mkdir -p $(dir $@)
	bwa index -p $(basename $@) ../$^

#Map reads against contigs
mapped.sam: $(R1) $(R2) | bwa_idx/iva_ctg7.bwt
	bwa mem -t $(threads) -M -T 40 -B 8 $(basename $|) $^ > $@

#Convert from sam to bam removing secondary mappings
%_sort.bam: %.sam
	samtools view -F 256 -hSb $^ | samtools sort -o $@ -T tmp -@ 8 -

#Convert from sam to bam removing secondary mappings
%_sort.bam: %.sam
	samtools view -F 256 -hSb $^ | samtools sort -o $@ -T tmp -@ 8 -

#Remove duplicates
%_sort_rmdup.bam : %_sort.bam
	run_picard MarkDuplicates INPUT=$^ OUTPUT=$@ METRICS_FILE=$*.metrics REMOVE_DUPLICATES=True ASSUME_SORTED=True

#Index resulting bamfile
%.bam.bai: %.bam
	samtools index $^

#***********************************************************
# Contig annotation with blast
#***********************************************************
nt_blastdb:=/labcommon/db/blastdb/nt/nt
nr_blastdb:=/labcommon/db/blastdb/nr/nr

%_tblastx_nt.txt: %.fa
	tblastx -evalue 0.001 -num_threads 14 -max_target_seqs 10 -outfmt 0 -show_gis -db $(nt_blastdb) -query $< -out $@

%_blastx_nr.txt: %.fa
	blastx -evalue 0.001 -num_threads 16 -max_target_seqs 10 -outfmt 0 -show_gis -db $(nr_blastdb) -query $< -out $@

#*************************************************************
# Variant calling with Freebase
#*************************************************************
%.vcf: %.bam ctgs/iva_ctg7.fa
	$(FREEBAYES_BIN) -f $(word 2,$^) -C 5 -p 1 --theta 0.05 \
		--min-mapping-quality 40 --min-base-quality 30 \
		 --pooled-continuous --report-monomorphic $< > $@
