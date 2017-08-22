#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $usage = "
remove_seqs.pl <id_list> <fasta file> <new output fasta file>

given a list of IDs in a file (one per line), this program will remove the
seqs with any of those IDs from the sequence file
";
die $usage unless @ARGV == 3;


open PAT, '<', $ARGV[0] or die "couldn't open $ARGV[0]: $!\n";
my @pats = <PAT>;
chomp @pats;
my %patterns = map { $_ => 1 } @pats;
close PAT;
undef @pats;

my $out = Bio::SeqIO->new('-file'   => ">$ARGV[2]",
			  '-format' => 'fasta');

my $db  = Bio::SeqIO->new('-file'   => $ARGV[1],
			  '-format' => 'fasta');


while (my $seq = $db->next_seq) {
    my $id = $seq->primary_id();
    if ( exists $patterns{$id} ) {
	print STDERR $id, "\n";
	next;
    }
    else {
	$out->write_seq($seq);
    }
}
