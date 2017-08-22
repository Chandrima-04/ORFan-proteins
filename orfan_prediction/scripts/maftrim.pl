#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2010-06-03

use Modern::Perl;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Bio::AlignIO;

my ($len);
GetOptions( 'len:i' => \$len );

my $usage = "
maftrim.pl - trim MAF aligments to a given length and write out in ClustalW

--len <n>   where n is an integer length

Alignments already shorter than the requested length will be skipped.

";
die $usage unless ($len && @ARGV);

foreach my $file (@ARGV) {

    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );

    my $in = Bio::AlignIO->new(
        '-file'   => $file,
        '-format' => 'maf'
    );
    my $outbase = $base . '_' . $len . '.aln';
    my $outfile = File::Spec->catfile( $path, $outbase );
    my $out     = Bio::AlignIO->new(
        '-file'   => ">$outfile",
        '-format' => 'clustalw',
        '-displayname_flat' => 1,
    );

    while ( my $aln = $in->next_aln() ) {
        if ( $aln->length >= $len ) {
            my $shortaln = $aln->slice( 1, $len );
            $out->write_aln($shortaln);
        }
    }
}
