#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-05-09

use Modern::Perl;
#use Getopt::Long;
# use File::Basename;
# use File::Spec;
# use Cluster;
# use ClusterGroup;
# use ClusterUtils;

# my ();
# GetOptions();

my $usage = "picklist2mcltab.pl - make mcl-style tabfile from a two-column picklist

Usage: picklist2mcltab.pl <picklist>

where
mcl tabfile is one line per cluster, tab-delimited seqids
picklist is family in the first column, seqid in the second.
e.g.
ATP_bind_3      EBC59051.1
ATP_bind_3      ECJ67731.1
ATP_bind_3      EBM82019.1
AAA_3   EBH79597.1
AAA_3   EDD97894.1

";
die $usage unless @ARGV;

# # parse filename
# my ( $base, $path, $suffix ) = fileparse( $ARGV[0], qr/\.[^.]*/ );
# my $outbase = $base . '.tab';
# my $outfile = File::Spec->catfile( $path, $outbase );

open(my $pick, '<', $ARGV[0]) or die "couldn't open $ARGV[0]:!$\n";

# put picklist into hash
my %families;
while (my $line = <$pick>) {
    chomp $line;
    my ($family, $seqid) = split /\t/, $line;
    push @{ $families{$family} }, $seqid;
}

# loop through hash, writing output as we go.
while (my ($famname, $seqids) = each %families) {
    print join("\t", @$seqids), "\n";
}