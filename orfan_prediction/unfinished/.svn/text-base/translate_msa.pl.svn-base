#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-05-12

use Modern::Perl;
use Getopt::Long;
use File::Basename;
use File::Spec;
use Bio::AlignIO;
use Bio::SeqIO;

my ($frame);
GetOptions('frame:i' => \$frame);

my $usage = "translate_msa.pl --frame <reading frame>

";
die $usage unless @ARGV && $frame;

foreach my $file (@ARGV) {

    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );

    my $in = Bio::AlignIO->new(
        '-file'   => $file,
        '-format' => 'clustalw',
    );

    # my $i = 0;
    # my $outbase = $base . '_' . $i++ . '.aln';
    # my $outfile = File::Spec->catfile( $path, $outbase );
    # my $out     = Bio::AlignIO->new(
    #     '-file'             => ">$outfile",
    #     '-format'           => 'fasta',
    #     '-displayname_flat' => 1,
    # );

    my $out = Bio::SeqIO->new('-fh'     => \*STDOUT,
                                '-format' => 'fasta');

    while ( my $aln = $in->next_aln() ) {

        foreach my $seq ( $aln->each_seq() ) {
            # flip negative frames first
            if ($frame =~ /^-(\d)/) {
                $seq = $seq->revcom();
                # my $sequence = $seq->seq;
                # $seq->seq(reverse($sequence));
                $frame = $1;
            }
            my $prot_seq = $seq->translate('-frame' => $frame);
            $out->write_seq($prot_seq);
        }

    }
}
