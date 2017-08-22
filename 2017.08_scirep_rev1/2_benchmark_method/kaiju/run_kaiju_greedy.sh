#!/bin/bash
set -euo pipefail

module load kaiju/1.5.0

KAIJU_DB_DIR=/labcommon/db/kaiju/progenomes_refseqvir
IN_FASTA=../../0_data/454_seqs/myunknown_goodclean5_prinsgood.fa

OUT_PREFIX=kaiju_greedy/454_seqs_kaiju_greedy

mkdir -p kaiju_greedy

kaiju -x -m 20 -a greedy -e 1 -z 16 -t ${KAIJU_DB_DIR}/nodes.dmp -f ${KAIJU_DB_DIR}/kaiju_db.fmi -i ${IN_FASTA} -o ${OUT_PREFIX}.txt

addTaxonNames -u -r superkingdom,family,species -t ${KAIJU_DB_DIR}/nodes.dmp -n ${KAIJU_DB_DIR}/names.dmp -i ${OUT_PREFIX}.txt -o ${OUT_PREFIX}.names.txt

kaijuReport -t ${KAIJU_DB_DIR}/nodes.dmp -n ${KAIJU_DB_DIR}/names.dmp -i ${OUT_PREFIX}.txt -r family -o ${OUT_PREFIX}.family.report
kaijuReport -t ${KAIJU_DB_DIR}/nodes.dmp -n ${KAIJU_DB_DIR}/names.dmp -i ${OUT_PREFIX}.txt -r genus -o ${OUT_PREFIX}.genus.report
kaijuReport -t ${KAIJU_DB_DIR}/nodes.dmp -n ${KAIJU_DB_DIR}/names.dmp -i ${OUT_PREFIX}.txt -r species -o ${OUT_PREFIX}.species.report

#Get virus hits
grep "Viruses" ${OUT_PREFIX}.names.txt > ${OUT_PREFIX}.names.virus.txt
