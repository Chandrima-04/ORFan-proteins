#!/usr/bin/perl
# $Id: evolver2clustal.pl 1916 2010-05-27 11:57:04Z dmessina $
# created by Dave Messina on 2010-05-27

use Modern::Perl;
use Bio::AlignIO;
use File::Basename;
use File::Spec;
use Getopt::Long;

my $outdir;
GetOptions( 'outdir:s' => \$outdir );

my $usage = "
evolver2clustal.pl - take multi-phylip files from evolver and
                     split them into single-clustalw files.
                     
Usage: evolver2clustal.pl --outdir <output dir>   evolver1.phy ... evolverN.phy

where

outdir is the directory where the output clustalw files should go.

Note: the output files will be named based on the input filename.

e.g. if input filename is cluster1.phy, output filenames will be
cluster1_1.aln, cluster1_2.aln, and so on.
";
die $usage                                            unless @ARGV;
die "\nno outdir!\n$usage"                            unless defined $outdir;
die "$outdir is not a directory or does not exist!\n" unless -d $outdir;

foreach my $file (@ARGV) {
    next unless $file =~ /\.phy$/;

    # make sure we can reach that file from anywhere on the filesystem
    $file = File::Spec->rel2abs($file);

    # grab file basename
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );

    # read evolver file (multi-phylip)
    my $in = Bio::AlignIO->new(
        '-file'        => $file,
        '-format'      => 'phylip',
        '-interleaved' => 0,
    );

    my $rec_number = 0;
    while ( my $aln = $in->next_aln() ) {
        $rec_number++;
        my $outbase = $base . '_' . $rec_number . '.aln';
        my $outfile = File::Spec->catfile( $outdir, $outbase );
        my $out     = Bio::AlignIO->new(
            '-file'   => ">$outfile",
            '-format' => 'clustalw',
            '-displayname_flat' => 1,
        );
        $out->write_aln($aln);
    }
}
