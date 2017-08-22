#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-04-04

use Modern::Perl;
use Getopt::Long;

#my ();
#GetOptions();

my $usage = "mcltab2id.pl - extract sequence IDs from an mcl tab file
";
die $usage unless @ARGV;

while (my $line = <>) {
    chomp $line;
    my @ids = split /\t/, $line;
    foreach my $member (@ids) {
        print STDOUT $member, "\n";
    }
}