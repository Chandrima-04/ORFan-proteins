#!/usr/bin/perl
# $Id: rose_fasta_split.pl 1676 2010-03-23 14:55:11Z dmessina $

use strict;
use warnings;
use Bio::SeqIO;

my $usage = "
rose_fasta_split.pl - split a multi-run ROSE fasta file into one per run

rose_fasta_split.pl <seqs_per_run> <rose.fa>
";
die $usage unless @ARGV == 2;

my ($seqcount, $infile) = @ARGV;

my $stream = Bio::SeqIO->new('-format' => 'fasta',
                             '-file'   =>  $infile,);

my $seqindex  = 0;
my $fileindex = 0;
my ($outname, $outfh);
while (my $seqobj = $stream->next_seq) {
    
    if ($seqindex % $seqcount == 0) {
        $fileindex++;
        $outname = $infile . ".$fileindex";
        $outfh = Bio::SeqIO->new('-format' => 'fasta',
                                 '-file'   => ">$outname");
    }

    $outfh->write_seq($seqobj);
    $seqindex++;
}

if ($seqindex % $seqcount != 0) {
    print STDERR "ERROR: last sequence index $seqindex not divisible by $seqcount\n";
}