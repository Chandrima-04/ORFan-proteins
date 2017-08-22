#!/usr/bin/perl
# $Id: parsed_paml_split.pl 1595 2010-02-02 08:20:51Z dmessina $

use warnings;
use strict;
use Test::Deep::NoTest;

my $usage = "
parsed_paml_split - compare duplicate measurements in a tab-delimited parsed
paml file and outputs a non-duplicated set.

lines should look like:
clustername    seq1    seq2    Ka/Ks    Ka    Ks
";

my %seen;

while (my $line = <>) {
    chomp $line;
    my ($cluster, $seq1, $seq2, $kaks, $ka, $ks) = split "\t", $line;
    die "bad line :$line" unless ($cluster && $seq1 && $seq2 && $kaks && $ka
        && $ks);
    
    # join the cluster and 2 seq names to make a unique hash key
    my $key = join('|', $cluster, $seq1, $seq2);

    # anon hash to store the other values
    my $value = { 'kaks' => $kaks,
                  'ka'   => $ka,
                  'ks'   => $ks,
                };
    
    # if it's a duplicate line, compare to stored values
    if ($seen{$key}) {
        
        my ($ok, $stack) = eq_deeply($seen{$key}, $value);
    
        # report details of any inconsistent data
        unless ($ok) {
            print STDERR "\n$stack\n";
        }
    }
    # if it's the first time we've seen data for this key, store it
    # and print it
    else {
        $seen{$key} = $value;
        print STDOUT $line, "\n";
    }
}