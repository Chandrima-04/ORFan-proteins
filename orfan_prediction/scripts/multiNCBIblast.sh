#!/bin/bash

module add blast || echo "ERROR: module call failed."

# vars
#     EXEC=blastp
#       DB=/afs/pdc.kth.se/home/d/dmessina/dave1/virus/test/test_orfs.fa
# QUERYDIR=/afs/pdc.kth.se/home/d/dmessina/dave1/virus/test
#    QUERY=test_orfs.fa
#  OPTIONS="E=1e-1 -cpus 8 nonnegok novalidctxok shortqueryok V=5000 B=5000 -wordmask seg+xnu"
#   OUTDIR=/afs/pdc.kth.se/home/d/dmessina/dave1/virus/test
#  OUTFILE=multiWUBlast_$$.out

# get vars from calling script
    EXEC=$5
      DB=$4
   DBDIR=$6
   QUERY=$1
QUERYDIR=$2
# OPTIONS="-m 8 -a 8"
 OPTIONS="-evalue 1e-4 -outfmt 7 -num_descriptions 20 -num_alignments 20 -db_soft_mask 30"
  OUTDIR=$3
 OUTFILE=$QUERY.out
  OUTERR=$QUERY.err

# set up working dir and copy files over
now=`date`
echo "start copying files to node - $now"
SCRATCH=/scratch/NCBIblasting
rm -rf $SCRATCH
mkdir $SCRATCH
cp $QUERYDIR/$QUERY $SCRATCH
cp $DBDIR/$DB.phd $SCRATCH
cp $DBDIR/$DB.phi $SCRATCH
cp $DBDIR/$DB.phr $SCRATCH
cp $DBDIR/$DB.pin $SCRATCH
cp $DBDIR/$DB.pnd $SCRATCH
cp $DBDIR/$DB.pni $SCRATCH
cp $DBDIR/$DB.ppd $SCRATCH
cp $DBDIR/$DB.ppi $SCRATCH
cp $DBDIR/$DB.psd $SCRATCH
cp $DBDIR/$DB.psi $SCRATCH
cp $DBDIR/$DB.psq $SCRATCH
now2=`date`
echo "done copying files to node - $now2"

# execute command
# ($EXEC -p blastp -d $SCRATCH/$DB -i $SCRATCH/$QUERY $OPTIONS ) 1> $SCRATCH/$OUTFILE 2> $SCRATCH/$OUTERR &
($EXEC -db $SCRATCH/$DB -query $SCRATCH/$QUERY $OPTIONS -out $SCRATCH/$OUTFILE ) 2> $SCRATCH/$OUTERR &
# give a little output
echo "query:$QUERY outfile:$OUTFILE"


# wait for jobs to finish
wait
now3=`date`
echo "done with jobs - $now3"

# copy back to my disk
cp $SCRATCH/$OUTFILE $OUTDIR
cp $SCRATCH/$OUTERR $OUTDIR
now4=`date`
echo "completely done - $now4"