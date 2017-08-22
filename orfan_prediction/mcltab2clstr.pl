#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-04-04

use Modern::Perl;
use Getopt::Long;

#my ();
#GetOptions();

my $usage = "mcltab2clstr.pl - convert mcl tab file to fake CD-HIT clstr format
";
die $usage unless @ARGV;

my $cluster_id = 0;
while (my $line = <>) {
    print STDOUT '>Cluster ', $cluster_id++, "\n";
    chomp $line;
    my @ids = split /\t/, $line;
    my $member_id = 0;
    foreach my $member (@ids) {
        print STDOUT $member_id++, "\t0aa, >", $member, ' at 100%', "\n";
    }
}