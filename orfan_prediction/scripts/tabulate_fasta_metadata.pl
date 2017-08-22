#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-08-17

use Modern::Perl;

my $usage = "tabulate_fasta_metadata.pl - extract details from metagenomic read headers
and write out a summary table

Usage: tabulate_fasta_metadata.pl 454.all.fasta
";
die $usage unless @ARGV;

my %samples;

while (my $line = <>) {
    next unless $line =~ /^>/;
    chomp $line;
    
    my @raw_fields = split /\s+/, $line;
    my %hash_fields;
    
    my $seq_id;
    if (substr($raw_fields[0], 0, 1) eq '>') {
        my $first = shift @raw_fields;
        $seq_id = substr $first, 1;
        
    }
    else {
        die "wha? first field is ", $raw_fields[0], "\n";
    }
    
    foreach my $field (@raw_fields) {
        my @pieces = split '=', $field;
        die "problem! $field\n" if scalar @pieces != 2;
        $hash_fields{$pieces[0]} = $pieces[1];
    }
    
    # origin
    my $origin = $hash_fields{'origin'} || 'none';
    
    # samples
    my $raw_sample = $hash_fields{'samples'} || 'none';
    my @samples = split ',', $raw_sample;

    foreach my $entry (@samples) {
        my $sample;
        if ( $entry =~ /(\S+)\(\d+\)/ ) {
            $sample = $1;
        }
        else { $sample = $entry; }
        
        $samples{$sample}->{'origin'}->{$origin}++;
    }


    # if (exists $samples{$sample}) { die "sample $sample already exists!\n"; }
    # else {
    #     $samples{$sample} = ('origin' => $origin,
    #                          # 'e'      => $evalue,
    #                          # 'nreads' => $nreads,
    #                          # 'type'   => \%type,
    #                         );
    # }
}

foreach my $entry (sort {$a cmp $b} keys %samples){
    my %origins = %{ $samples{$entry}->{'origin'} };
    foreach my $origin (sort keys %origins) {
        say join "\t", $entry, $origin, $origins{$origin};
    }
}

1;