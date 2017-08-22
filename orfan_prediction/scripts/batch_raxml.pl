#!/usr/bin/env perl
# $Id: batch_raxml.pl 1694 2010-03-31 18:43:02Z dmessina $

use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $usage = "
batch_raxml.pl - run RAxML on multiple alignments in one go

Usage: batch_raxmpl.pl cluster1.phylip ... clusterN.phylip

NOTE: This program assumes input alignments have a cluster[0-9]+ in their name
      which will be passed to RAxML for naming the output files.
      
";
die $usage unless @ARGV > 0;

foreach my $file (@ARGV) {

    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my $seqfile = $file;

    my $outbase;
    if ( $seqfile =~ /(cluster\d+)/ ) {
        $outbase = $1;
    }
    else {
        die "couldn't find cluster[0-9] regex in $seqfile filename!\n";
    }

    my $command =
"raxml -T 2 -f a -x 59457 -p 59457 -# 100 -m GTRGAMMA -s $seqfile -n $outbase";
    my $retval = system($command);

    if ( $retval > 0 ) { print STDERR $retval, "\t$!\n"; }
}
