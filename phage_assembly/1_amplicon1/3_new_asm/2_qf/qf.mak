# Raw reads quality filtering pipeline for MAARS samples

# Author: Mauricio Barrientos-Somarribas
# Email:  mauricio.barrientos@ki.se

# Copyright 2014 Mauricio Barrientos-Somarribas
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#Make parameters
SHELL := /bin/bash

sample_name := d9539_amplicon1

#Temporary merged R1 and R2
R1_TMP :=  $(wildcard ../reads/*_1.fastq.gz)
R2_TMP := $(wildcard ../reads/*_2.fastq.gz)

#Outfile
OUT_PREFIX := $(sample_name)_qf

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Uncomment for debugging, otherwise make deletes all intermediary files
.SECONDARY:

.PHONY: all

all: 4_phix_rm

#*************************************************************************
# Call to NESONI for quality trimming and clipping illumina adapters
#*************************************************************************
#You have to specify quality is phred33 because with cutadapt clipped fragments nesoni fails to detect encoding
1_nesoni: $(R1_TMP) $(R2_TMP)
	mkdir -p $@
	nesoni clip --adaptor-clip no --homopolymers yes --qoffset 33 --quality 5 --length 50 \
		--out-separate yes $@/$(OUT_PREFIX) pairs: $^

2_cutadapt: 1_nesoni
	mkdir -p $@
	#Remove Illumina TruSeq Barcoded Adapter from fwd and rev pair
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
			 -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
				--overlap=5 --error-rate=0.1 --minimum-length 50 \
				-o $@/$(OUT_PREFIX)_R1.fq.gz -p $@/$(OUT_PREFIX)_R2.fq.gz $</*_{R1,R2}.fq.gz
	#Remove reverse complement of Illumina TruSeq Universal Adapter from reverse pair
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
			 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
				--overlap=5 --error-rate=0.1 --minimum-length 50 -o $@/$(OUT_PREFIX)_single.fq.gz $</*_single.fq.gz

#*************************************************************************
# QC for nesoni
#*************************************************************************
2_cutadapt_qc: 2_cutadapt
	mkdir -p $@
	fastqc --noextract -k 10 -o $@ $^/*.fq.gz
	bash jellyfish2.sh $@/$(OUT_PREFIX) $^

#*************************************************************************
# Pair Merging with flash
#*************************************************************************
3_mergepairs: 2_cutadapt
	mkdir -p $@
	flash -m 10 -M 300 $^/*R1.fq.gz $^/*R2.fq.gz -o $(OUT_PREFIX) -d $@

#*************************************************************************
# Remove phix contamination
#*************************************************************************
4_phix_rm: 3_mergepairs
	mkdir -p $@/phix_bwa
	bwa index -p $@/phix_bwa/phix174 /labcommon/db/fasta/contaminants/phiX174.fasta
	bwa mem -t 16 -k 15 -M -T 30 $@/phix_bwa/phix174 $</*_1.fastq $</*_2.fastq \
	 	| samtools view -hSb -f 12 - | samtools bam2fq - | gzip > $@/$(OUT_PREFIX)_pe.fq.gz
	bwa mem -t 16 -k 15 -M -T 30 $@/phix_bwa/phix174 <(gzip -c 2_cutadapt/$(OUT_PREFIX)_single.fq.gz | cat $</*.extendedFrags.fastq) \
	 	| samtools view -hSb -f 4 - | samtools bam2fq - | gzip > $@/$(OUT_PREFIX)_se.fq.gz


.PHONY: clean
clean:
	@echo "Nothing to clean??"
