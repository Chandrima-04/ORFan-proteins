#!/usr/bin/perl
# $Id: kalign2clustalw.pl 1891 2010-05-24 13:28:30Z dmessina $
# created by Dave Messina on 2010-05-20

use Modern::Perl;
use File::Basename;

my $header = "CLUSTAL W (1.82) multiple sequence alignment\n";
my $outdir = 'cluster_aln2';

foreach my $file (@ARGV) {
    open(my $infh, '<', $file) or die "couldn't open $file: $!\n";
    
    my $base = basename($file);
    open(my $outfh, '>', "$outdir/$base") or die "couldn't open $outdir/$base: $!\n";
    print $outfh $header;
    
    while (<$infh>) {
        next if $. == 1;
        print $outfh $_;
    }
    close $outfh or die "couldn't close $outdir/$base: $!\n";
    close $infh  or die "couldn't close $file: $!\n";
}