#!/usr/bin/perl
# $Id: strip_nse.pl 3068 2011-07-06 11:51:37Z dmessina $
# created by Dave Messina on 2011-06-30

use Modern::Perl;
use Getopt::Long;
use Bio::AlignIO;
use File::Basename;
use File::Spec;

my ($indir, $outdir);
GetOptions('indir:s'  => \$indir,
           'outdir:s' => \$outdir);

my $usage = "strip_nse.pl - remove nse from clustalw alignments

strip_nse.pl <.aln files> --outdir <dir to write new .aln files>
";
die $usage unless @ARGV && -d $outdir;

foreach my $file (@ARGV) {

    my $aln_in = Bio::AlignIO->new(
        '-file'  => $file,
        '-format' => 'clustalw');

    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my $outbase = $base . $suffix;
    my $outfile = File::Spec->catfile($outdir, $outbase);
    my $aln_out = Bio::AlignIO->new('-file'  => ">$outfile",
                                    '-format' => 'clustalw',
                                    '-displayname_flat' => 1);

    while (my $aln = $aln_in->next_aln) {
        $aln_out->write_aln($aln);
    }
}
