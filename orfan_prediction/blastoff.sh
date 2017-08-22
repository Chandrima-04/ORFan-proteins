#!/bin/sh
# $Id: blastoff.sh 2060 2010-06-15 13:24:33Z dmessina $

# USAGE -- do this on ferlin headnode
# a07c01n08$ for f in {0..65}; do echo $f; esubmit -v -n 1 -t 1200 blastoff.sh $f; done
#   where 0..65 are the indexes of the shattered fasta files made by dbshatter

if [[ $# != 1 ]]; then
    echo "You must supply an argument!"
    exit 1
fi

QUERYDIR=~/dave1/virus/data/cluster_dna
OUTERR="cluster$1.err"

cd /scratch
now5=`date`
echo "start copying tarfiles - $now5" >> $OUTERR
rcp a08c31n10.pdc.kth.se:/scratch/uniref50_blastdb.tgz .
now6=`date`
echo "done copying tarfiles - $now6"
echo "done copying tarfiles - $now6" >> $OUTERR

now=`date`
echo "start untaring files - $now" >> $OUTERR
tar zxf uniref50_blastdb.tgz
echo "untaring done" >> $OUTERR

QUERY="$QUERYDIR/cluster$1.fa"
cp $QUERY .

EXEC=~/bin/ncbi-blast-2.2.23+_x64/bin/blastx
OPTIONS="-evalue 1e-4 -outfmt 7 -num_descriptions 20 -num_alignments 20 -num_threads 8"
OUTFILE="cluster$1.uniref50"

RESULTDIR=~/dave1/virus/20100615/vs_uniref50
now2=`date`
echo "starting blast - $now2" >> $OUTERR
($EXEC -db uniref50 -query $QUERY $OPTIONS -out $OUTFILE ) 2>> $OUTERR &
now3=`date`
echo "finished blast - $now3" >> $OUTERR

wait
cp $OUTFILE $RESULTDIR
now=`date`
echo "done copying to $OUTFILE to $RESULTDIR - $now4" >> $OUTERR
cp $OUTERR $RESULTDIR
