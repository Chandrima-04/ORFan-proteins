#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $in = Bio::SeqIO->new('-format' => 'fasta',
			 '-file'   => $ARGV[0],);

my $out = Bio::SeqIO->new('-format' => 'fasta',
			  '-fh'     => \*STDOUT,);

while (my $seqobj = $in->next_seq) {
    my $seq = $seqobj->seq;
    my $len = length($seq);
    my $lc_count = $seq =~ tr/a-z/a-z/;
    if ($len-$lc_count < 80) {
	next;
    }
    else {
	$out->write_seq($seqobj);
    }
}
