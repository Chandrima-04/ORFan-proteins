SHELL:=/bin/bash

sample_name := d9539_amplicon1

paired_end:=../2_qf/4_phix_rm/$(sample_name)_qf_pe.fq.gz
single:=../2_qf/4_phix_rm/$(sample_name)_qf_se.fq.gz

MASURCA_BIN:=/labcommon/tools/MaSuRCA-2.3.2/bin

ifndef threads
	threads := 16
endif

#Delete produced files if step fails
.DELETE_ON_ERROR:

#Avoids the deletion of files because of gnu make behavior with implicit rules
.SECONDARY:

all: fermi megahit spades

#*****************************************************************************
# ASSEMBLY
#*****************************************************************************
#Megahit
megahit: $(paired_end) $(single)
	megahit -m 5e10 -l 600 --k-step 2 --k-max 101 --input-cmd "cat $^" --cpu-only -t $(threads) -o $@

#Spades
spades: $(paired_end) $(single) $(single)
	mkdir -p $(dir $@)
	spades.py -t $(threads) --12 $< -s $(word 2,$^) -o $@ -m 64 -k $$(seq -s , 21 4 101)


fermi: $(paired_end) $(single)
	mkdir -p $@
	cd $@ && run-fermi.pl -t $(threads) -k 40 $(addprefix ../,$^) > fermi.mak
	cd $@ && make -f fermi.mak -j $(threads) fmdef.p2.mag.gz
	seqtk seq -A $@/fmdef.p2.mag.gz > $@/fermi_contigs.fa


#Sanger IVA with all qf reads
iva_asm: $(paired_end) $(single)
	source activate py3k && iva -f <(seqtk seq -1 $<) -r <(seqtk seq -2 $<) --contigs $(word 2,$^) -t $(threads) iva_asm

#*************************************************************
# Other assemblers
#*************************************************************
pe_reads/masurca_R%.fq.gz: $(paired_end)
	mkdir -p $(dir $@)
	seqtk seq -$* $< | gzip > $@

masurca: masurca-qf.config pe_reads/masurca_R1.fq.gz pe_reads/masurca_R2.fq.gz
	$(MASURCA_BIN)/masurca $<
	bash assemble.sh

#**************
# Clean
#**************

.PHONY: clean

clean:
	-rm -r pe_reads
