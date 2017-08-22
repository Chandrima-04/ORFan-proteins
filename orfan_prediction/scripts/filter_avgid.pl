#!/usr/bin/perl
# $Id: filter_avgid.pl 3108 2011-09-15 07:56:55Z dmessina $
# created by Dave Messina on 2011-09-14

use Modern::Perl;
use File::Basename;

my $avgid_file = shift @ARGV;
open (my $avgid, "<", $avgid_file) or die "couldn’t open $avgid_file: $!";

# read in avg id info
while (my $line = <$avgid>) {
    chomp $line;
    my ($file, $label, $rawvalue) = split /:/, $line;
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    
    # this is the actual avg % id number
    my $value;
    if ($rawvalue =~ /\s+(\d+)\%/) {
        $value = $1;
    }
    
    # filtering out anything 92% id or higher
    next if $value >= 92;
    
    say $base;
}
