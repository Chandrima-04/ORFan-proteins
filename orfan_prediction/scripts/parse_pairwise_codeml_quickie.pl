#!/usr/bin/perl
# $Id: parse_pairwise_codeml_quickie.pl 1827 2010-05-06 09:08:57Z dmessina $

use strict;
use warnings;
use File::Basename;

my $usage = "
parse_pairwise_codeml_quickie - grabs Ka/Ks info from codeml result file(s)
                                and reports in a tab-delimited format.

Report format is:
Cluster    Seq1    Seq2    dN/dS    dN    dS

NOTE: This program assumes the codeml result files are in files or folders
named after their cluster, like 'cluster123'.

";
die $usage unless @ARGV;

while (my $infile = shift @ARGV) {
    my ( $name, $path, $suffix ) = fileparse( $infile, qr/\.[^.]*/ );

    # extract cluster name from infile (including path)
    my ($cluster) = $infile =~ /(cluster\d+)/;
    die "no cluster name in $infile\n" unless $cluster;

    # slurp file
    open(my $fh, $infile) or die "couldn't open $infile";
    my @lines = <$fh>;
    chomp @lines;

    # success flag
    my $success;
    
    # grab the line which looks like this:
    # all_c13442/1-153     0.8448 (1.1901 1.4087)
    LINE: for (my $i = $#lines; $i >= 0;  $i--) {
        my $line = $lines[$i];

       # validate it
       if ($line =~ m/
                             (.+?)              # name of second seq in pair
                             \s+                # spaces after sequence name
                             (-?(?:\d+\.)?\d+)  # dNdS
                             \s                 # space
                             \(                 # opening paren
                             (-?(?:\d+\.)?\d+)  # dN
                             \s                 # space
                             (-?(?:\d+\.)?\d+)  # dS
                             \)                 # closing paren
                             $                  # end of line
                             /x ) {
    
            my ($seq2, $omega, $dN, $dS) = ($1, $2, $3, $4);
    
            # grab seq1's name off previous line
            my $seq1_line = $lines[$i-1];
            my ($seq1) = $seq1_line =~ /(.+?)\s+$/;
    
            # report results if we got everything
            if ($seq1 && $seq2 && $omega && $dN && $dS) {
                print STDOUT join("\t", $cluster, $seq1, $seq2, $omega, $dN, $dS), "\n";
                $success = 1;
                last LINE;
            }
        }
    }

    unless ($success) {
        print STDERR "couldn't find dNdS in file $infile\n";
    }
}
