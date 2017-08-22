#!/bin/bash


for aln in $@; do
    g=`basename $aln ".aln"`;
    echo $g;
    outfile=${g}_rnacode.txt
    RNAcode --outfile $outfile --eps --eps-cutoff 1.0 --eps-dir ${g}_eps --tabular $aln
done
