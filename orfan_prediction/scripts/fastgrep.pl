#!/usr/bin/perl

use strict;
use warnings;

my $usage = "
fastgrep.pl <id_list> <tab blast file> <new output file>

given a list of IDs in a file (one per line), this program will look for those IDs to be
queries or hits in the tabbed blast file (in the first or second columns).

Any lines with one of the IDs in the list will be SKIPPED.

";
die $usage unless @ARGV == 3;

open PAT, '<', $ARGV[0] or die "couldn't open $ARGV[0]: $!\n";
my @pats = <PAT>;
chomp @pats;
my %patterns = map { $_ => 1 } @pats;
close PAT;

open OUT, '>', $ARGV[2] or die "couldn't open $ARGV[2]: $!\n";

open DB, '<', $ARGV[1] or die "couldn't open $ARGV[1]: $!\n";
while (my $line = <DB>) {
    if ($line =~ /^#/) {
	print OUT $line;
	next;
    }

    my ($query, $hit) = split /\t/, $line, 3;

    next if ( exists($patterns{$query}) );
    next if ( exists($patterns{$hit})   );

    print OUT $line;
}
