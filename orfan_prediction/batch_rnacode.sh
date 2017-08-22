#!/bin/bash

outroot=$2

for aln in `ls $1`
do
    g=`basename $aln ".aln"`;
    echo $g;
    outfile=$outroot/${g}_rnacode.txt
    echo $outfile;
    RNAcode --outfile $outfile --eps --eps-cutoff 1.0 --eps-dir $outroot/${g}_eps --tabular $1/$aln
done
