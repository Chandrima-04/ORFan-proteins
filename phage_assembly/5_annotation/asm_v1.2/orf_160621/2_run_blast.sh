#!/bin/bash
samplename=d9539_asm_v1.2_orf

echo "Processing ${samplename}"
make -r -f 2_blast.mak sample_name=${samplename} nr
make -r -f 2_blast.mak sample_name=${samplename} env_nr
make -r -f 2_blast.mak sample_name=${samplename} hmp_pep
make -r -f 2_blast.mak sample_name=${samplename} metahit_pep
make -r -f 2_blast.mak sample_name=${samplename} metahit_2014_pep
echo
