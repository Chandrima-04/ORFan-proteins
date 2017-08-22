#!/usr/bin/perl

use strict;
use warnings;
my @ten;
my @five;
my @two;

while (my $line = <>) {
    chomp $line;
    my ($file, $size) = split ' ', $line;
    if ($size >= 10) {
	push @ten, $file;
    }
    elsif ($size >=5) {
	push @five, $file;
    }
    elsif ($size >= 2) {
	push @two, $file;
    }
    else {
	print STDERR  "whoa! small fry $file,\n";
    }
}

open my $tenfh, ">ten.txt" or die "couldn't ten:!$\n";
print $tenfh join("\n", @ten), "\n";
close $tenfh;

open my $fivefh, ">five.txt" or die "couldn't five:!$\n";
print $fivefh join("\n", @five), "\n";
close $fivefh;

open my $twofh, ">two.txt" or die "couldn't two:!$\n";
print $twofh join("\n", @two), "\n";
close $twofh;
