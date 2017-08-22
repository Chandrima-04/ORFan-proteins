#!/bin/bash
set -euo pipefail

module load samtools/1.3

#PICARD_JAR=/labcommon/tools/picard-tools/2.1.1/picard.jar

IN_BAM=../1_map_to_phage/454_seqs_to_phage.bam

samtools fasta ${IN_BAM} > 454_seqs_phage.fa
#java -jar ${PICARD_JAR} SamToFastq IN=${IN_BAM} OUT=
