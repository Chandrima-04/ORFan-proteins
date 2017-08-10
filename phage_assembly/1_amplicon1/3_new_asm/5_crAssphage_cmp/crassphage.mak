SHELL:=/bin/bash

sample_name:= d9539_amplicon1

threads := 16

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: $(addsuffix .bam.bai,ctgs/megahit ctgs/fermi )

#*********************************************************
# Create a last indexes from the crAssphage reference
#*********************************************************
crassphage/last_idx/crassphage.prj: crassphage/crAssphage.fasta
	mkdir -p $(dir $@)
	lastdb -cR01 $(basename $@) $^

#*************************************************************
# Align amplicon2 contigs against amplicon1 contig 7
#*************************************************************
ctgs/megahit.maf: ../3_new_asm/3_asm/megahit/final.contigs.fa | crassphage/last_idx/crassphage.prj
ctgs/fermi.maf:   ../3_new_asm/3_asm/fermi/fermi_contigs.fa   | crassphage/last_idx/crassphage.prj
# ctgs/spades.maf:  ../3_new_asm/spades/contigs.fasta     | amplicon1_iva_ctg7/last_idx/iva_ctg7.prj

%.maf:
	mkdir -p $(dir $@)
	parallel-fasta -j 8 "lastal $(basename $| )" < <(seqtk seq -L 500 $< ) > $@

%.bam: %.maf | sam.header
	maf-convert sam $< | cat $| - | samtools view -hSb - | samtools sort -O bam -T sort_tmp -@ $(threads) -o $@ -

%.bam.bai : %.bam
	samtools index $^
