#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-08-22

use Modern::Perl;
use Bio::SeqIO;
use File::Basename;
use File::Spec;

my $usage = "X2dash.pl - convert Xs to dashes '-' so that Epipe won't choke
Usage: X2dash.pl <aligned.fasta>

Note! should be run from the dir where you want output files to be written
";

foreach my $file (@ARGV) {
    my $seqio = Bio::SeqIO->new('-format' => 'fasta',
                                '-file'   => $file,
                                '-alphabet' => 'protein');

    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my $outfile = $base . '_noX.fa';
    my $seqout = Bio::SeqIO->new('-format' => 'fasta',
                                 '-file'   => ">$outfile",
                                 '-alphabet' => 'protein');

    while (my $seqobj = $seqio->next_seq) {
        my $seq = $seqobj->seq();

        say STDERR $seqobj->display_id(), "\t", length($seq);

        die ("too short before! ", $seqobj->display_id()) if $seqobj->length == 0;
        $seq =~ s/X/-/g;
        $seq =~ s/\*/-/g;
        die ("too short after! ", $seqobj->display_id()) if $seqobj->length == 0;
        $seqobj->seq($seq);
        $seqout->write_seq($seqobj);
    }
}