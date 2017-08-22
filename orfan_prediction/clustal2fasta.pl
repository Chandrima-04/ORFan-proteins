#!/usr/bin/perl

use strict;
use warnings;
use Bio::AlignIO;
use File::Basename;

my $usage = 'clustal2fasta - convert multiple alignments in clustalw format to
             multiple alignments in fasta format

usage: clustal2fasta <foo.aln> [bar.aln]

output file foo.fasta will be automatically created
';
@ARGV or die $usage;

foreach my $infile (@ARGV) {
    my $informat = 'clustalw';
    my $in       = Bio::AlignIO->new(
        -file   => $infile,
        -format => $informat,
    );

    # make output filename based on input filename
    my ( $base, $path, $suffix ) = fileparse( $infile, qr/\.[^.]*/ );
    my $outformat = 'fasta';
    my $outfile   = $path . $base . qq{.} . $outformat;

    my $out = Bio::AlignIO->new(
        -file   => ">$outfile",
        -format => $outformat,
        -interleaved => 0,
        -idlength    => 256,
    );

	while (my $aln = $in->next_aln()) {
		$out->write_aln($aln);
	}
}
