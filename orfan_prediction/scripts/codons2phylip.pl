#!/usr/bin/perl
# $Id: codons2phylip.pl 1678 2010-03-25 15:04:41Z dmessina $

use strict;
use warnings;
use Bio::AlignIO;
use File::Basename;

my $usage = "codons2phylip - convert sorta-phylip format from PAML's evolver to
                normal unbroken phylip format

usage: codons2phylip <foo.phylip> [bar.phylip]

output file foo.phy will be automatically created
";
@ARGV or die $usage;

foreach my $infile (@ARGV) {
    my $informat = 'phylip';
    my $in       = Bio::AlignIO->new(
        -file   => $infile,
        -format => $informat,
        -interleaved => 0,
    );

    # make output filename based on input filename
    my ( $base, $path, $suffix ) = fileparse( $infile, qr/\.[^.]*/ );
    my $outformat = 'phylip';
    my $outfile   = $path . $base . '_new.phy';

    my $out = Bio::AlignIO->new(
        -file        => ">$outfile",
        -format      => $outformat,
        -interleaved => 0,
        -idlength    => 30,
    );

    while ( my $aln = $in->next_aln() ) {
        $out->write_aln($aln);
    }
}