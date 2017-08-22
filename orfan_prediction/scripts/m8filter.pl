#!/usr/bin/perl
# $Id: m8filter.pl 1608 2010-02-12 16:02:25Z dmessina $

use strict;
use warnings;
use Getopt::Long;

my ($cutoff);
GetOptions('cutoff:s' => \$cutoff);

my $usage =
"m8filter.pl - filter an NCBI-BLAST m8 (tab-delimited) file on an e-value cutoff
               
Usage: m8filter.pl --cutoff 2.5e-05 myblast.m8

";
die $usage unless @ARGV && $cutoff;

while (my $line = <> ) {
    chomp $line;
    
    # m8 format fields are:
    # query_id, subject_id,
    # %ID, alignment length, # of mismatches, # of gaps, query_start,
    # query_end, subject_start, subject_end, e-value, bit score
    my @fields = split(/\t/, $line);
    
	print $line, "\n" if $fields[10] < $cutoff
}
