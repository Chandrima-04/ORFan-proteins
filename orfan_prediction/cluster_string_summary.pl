#!/usr/bin/perl
# $Id: cluster_string_summary.pl 2942 2011-04-20 14:05:13Z dmessina $
# created by Dave Messina on 2011-04-20

use strict;
use warnings;

my $usage = '
cluster_string_summary.pl - summarize results of clustering comparison, by cluster

$ cat *.bccout | cluster_string2stats.pl > bccout.stats
$ cluster_string_summary.pl bccout.stats
';
die $usage unless @ARGV;

my %clusters_from;
my %clusters_to;
my %clusters_bal;
while (my $line = <>) {
    chomp $line;
    next unless $line =~ /^MOVED/;
    my (undef, $seqid, $from, $to) = split "\t", $line;
    $clusters_from{$from}++;
    $clusters_to{$to}++;
}

foreach my $cluster ( sort { $a <=> $b } keys %clusters_from ) {
    my $to = defined $clusters_to{$cluster} ? $clusters_to{$cluster} : 0;
    $clusters_bal{$cluster} = abs($clusters_from{$cluster} - $to);
}

foreach my $cluster ( sort { $clusters_bal{$b} <=> $clusters_bal{$a} } keys %clusters_bal ) {

    my $to = defined $clusters_to{$cluster} ? $clusters_to{$cluster} : 0;
    print STDOUT join("\t", $cluster, $clusters_bal{$cluster},
        $clusters_from{$cluster}, $to), "\n";
}