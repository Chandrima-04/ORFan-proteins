#!/usr/bin/perl
# $Id: clustal_dumpfasta.pl 3070 2011-07-07 08:03:55Z dmessina $
# created by Dave Messina on 2011-07-06

use Modern::Perl;
use Bio::AlignIO;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;
use File::Spec;

my $outdir;
GetOptions('outdir:s' => \$outdir);

my $usage = "clustal_dumpfasta.pl - dump out sequences in a clustal MSA to
gapless unaligned FASTA seqs

Usage: clustal_dumpfasta.pl --outdir <outdir> *.aln

";
die $usage unless @ARGV && -d $outdir;

foreach my $file (@ARGV) {
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my $outbase = $base . '.fa';
    my $outfile = File::Spec->catfile($outdir, $outbase);

    my $aln_in = Bio::AlignIO->new(-format => 'clustalw',
                                   -file   => $file,
                                   -alphabet => 'protein',
                                   -displayname_flat => 1,
                                   -verbose => 1 );

    my $out = Bio::SeqIO->new(-format => 'fasta',
                              -file   => ">$outfile");

    my $aln = $aln_in->next_aln();
                          
    foreach my $seq ($aln->each_seq) {
        my $seq_string = $seq->seq;
        $seq_string =~ s/-//g;

        $seq->seq($seq_string);
        
        $out->write_seq($seq);
    }
}