#!/bin/bash
set -euo pipefail -o verbose
IFS=$'\n\t'

kraken_db=/labcommon/db/kraken/minikraken_20141208/
in_fasta=../../0_data/454_seqs/myunknown_goodclean5_prinsgood.fa

OUT_PREFIX=454_seqs_kraken
module load kraken/git_march_2016

kraken --preload --threads 16 --db ${kraken_db} --fasta-input --output ${OUT_PREFIX}.out ${in_fasta} 2> kraken.log

kraken-filter --db ${kraken_db} --threshold 0.05 ${OUT_PREFIX}.out > ${OUT_PREFIX}_filt.out

kraken-report --db ${kraken_db} ${OUT_PREFIX}_filt.out > ${OUT_PREFIX}_filt.report

#Subset virus hits from the report
cat 454_seqs_kraken_filt.report | awk 'BEGIN{flag=0} /Viruses/{flag=1} {if(flag==1)print $0}' > ${OUT_PREFIX}_filt.virus.report
