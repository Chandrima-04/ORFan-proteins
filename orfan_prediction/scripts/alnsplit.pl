#!/usr/bin/perl
# $Id: alnsplit.pl 1999 2010-06-09 09:22:28Z dmessina $
# created by Dave Messina on 2010-06-09

use Modern::Perl;
use Bio::AlignIO;
use Getopt::Long;
use File::Basename;
use File::Spec;

my $format = 'clustalw';
GetOptions( 'format:i' => \$format );

my $usage = "alnsplit.pl --format <clustalw> file1 ... fileN

--format   specify the format of the multiple alignment file you're splitting.
           default: clustalw
";
die $usage unless @ARGV;

foreach my $file (@ARGV) {

    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );

    my $in = Bio::AlignIO->new(
        '-file'   => $file,
        '-format' => $format
    );

    my $i = 0;
    while ( my $aln = $in->next_aln() ) {
        my $outbase = $base . '_' . $i++ . '.aln';
        my $outfile = File::Spec->catfile( $path, $outbase );
        my $out     = Bio::AlignIO->new(
            '-file'             => ">$outfile",
            '-format'           => $format,
            '-displayname_flat' => 1,
        );

        $out->write_aln($aln);
    }
}
