SHELL:=/bin/bash

R1:=$(wildcard reads/*_1.fastq.gz)
R2:=$(wildcard reads/*_2.fastq.gz)

SAMPLE:= sonnhamer

#Logging
log_name := $(CURDIR)pipe_$(shell date +%s).log
log_file := >( tee -a $(log_name) >&2 )

threads := 16

MASURCA_BIN:=/labcommon/tools/MaSuRCA-2.2.1/bin
#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: contigs.fa

#*************************************************************
# Quality trimming with nesoni
#*************************************************************
#Run quality filtering on the reads
nesoni/$(SAMPLE)_R1.fq.gz: $(R1) $(R2)
	mkdir -p nesoni
	nesoni clip --homopolymers yes --quality 20 --length 75 --out-separate yes \
			$(dir $@)/$(SAMPLE) pairs: $^ 2>> $(log_file)


#*************************************************************
# Subset reads that map to viral for assembly
#*************************************************************
#Create Refseq viral bwa index
bwa_idx/refseq_viral.% : /labcommon/db/fasta/viral.1.1.genomic.fna
	mkdir -p bwa_idx
	cd bwa_idx && bwa index -p refseq_viral $^

#Run mapping against refseq sequences
viral_subset/viral.sam: nesoni/$(SAMPLE)_R1.fq.gz nesoni/$(SAMPLE)_R2.fq.gz | bwa_idx/refseq_viral.bwt
	bwa mem -t $(threads) -M bwa_idx/refseq_viral $(wordlist 1,2,$^) > $@ 2>> $(log_file)

#Convert from sam to bam removing secondary mappings
%.bam: %.sam
	samtools view -F 256 -hSb -o $@ $^ 2> $(log_file)

#Sort files by coordinates
%_sortcoord.bam : %.sam
	run_picard SortSam.jar INPUT=$^ OUTPUT=$@ SORT_ORDER=coordinate 2> $(log_file)

#Remove duplicates
%_rmdup.bam: %.bam
	run_picard MarkDuplicates.jar INPUT=$^ OUTPUT=$@ METRICS_FILE=metrics.rmdup REMOVE_DUPLICATES=True

#Sort files by queryname
%_sortqn.bam : %.bam
	run_picard SortSam.jar INPUT=$^ OUTPUT=$@ SORT_ORDER=queryname 2> $(log_file)

#Keep only reads that mapped confidently (with both pairs)
%_mapped.bam: %.bam
	samtools view -f 2 -b $^ > $@ 2> $(log_file)

#Convert unmapped reads to Fastq for assembly
%_R1.fq: %.bam
	run_picard SamToFastq.jar INPUT=$^ FASTQ=$*_R1.fq SECOND_END_FASTQ=$*_R2.fq 2> $(log_file)

#*************************************************************
# Assembly
#*************************************************************
#Assemblies
.PHONY: abyss-qf masurca-qf masurca-virus iva-qf

#Masurca with only reads premapped to refseq viral
masurca_config.txt:
	$(MASURCA_BIN)/masurca -g $@

contigs.fa: masurca_config.txt viral_subset/viral_sortcoord_rmdup_sortqn_mapped_R1.fq
	$(MASURCA_BIN)/masurca $<
	bash assemble.sh

#Abyss with all qf reads
abyss-qf: ../nesoni/sonnhamer_R1.fq ../nesoni/sonnhamer_R2.fq
	mkdir -p abyss
	abyss-pe k=120 name=abyss np=16 in="$^"

#Sanger IVA with all qf reads
iva-raw: $(R1) $(R2)
	source activate py3k && iva -f $< -r $(word 2,$^) $@

