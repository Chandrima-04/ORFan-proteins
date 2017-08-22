#!/usr/bin/perl
# $Id: kaks_stats.pl 1833 2010-05-06 12:58:55Z dmessina $

use strict;
use warnings;
use feature 'switch';
use Statistics::Descriptive;
use Getopt::Long;

my ( $sim, $sort, $tab );
GetOptions(
    'sim'    => \$sim,
    'sort:s' => \$sort,
    'tab'    => \$tab,
);

my $usage = "
kaks_stats.pl - calculate distribution statistics on pairwise cluster Ka/Ks
                data
                
Usage: kaks_stats.pl cluster_kaks.txt

where cluster_kaks.txt is in tab-delimited format as produced by 
    parse_pairwise_codeml_quickie.pl
    
Clusters whose median < 0.5 have a '<' at the end of the line.
Clusters whose median ~= 1  have a '=' at the end of the line.
Clusters whose median > 1.5 have a '>' at the end of the line.

Options:
--sim   cluster names take the form cluster16fake.100 and should be combined
--sort  sort results by: 'cluster', 'median', 'stdev' (=default)
--tab   produce output in tab-delimited format (default is space-delimited)
";
die $usage unless @ARGV;

my ( %data, %errors );

# read data and group it by cluster
while ( my $line = <> ) {
    chomp($line);

    # skip comment lines
    next if $line =~ /^#/;

    # crude format check
    my @fields = split( /\t/, $line );
    my $field_count = scalar @fields;
    if ( $fields[1] eq '-no sites-' ) {
        print STDERR "WARNING: $fields[0] has no sites.\n";
        next;
    }
    elsif ( $field_count != 6 ) {
        die(
            "format on line $. doesn't look right; ",
            "expecting 6 fields, got $field_count.\n",
        );
    }

    # skip header
    if ( $. == 1 && $fields[3] eq 'dNdS' ) { next; }

    my ( $cluster, $seq1, $seq2, $dNdS, $dN, $dS ) = @fields;

    # extract cluster name from long name, e.g.
    # cluster16fake.100 => cluster16
    if ($sim) {
        my $full_cluster_name = $cluster;
        if ( $full_cluster_name =~ /(cluster\d+)fake\.(\d+)/ ) {
            $cluster = $1;
        }
    }

    # record and skip erroneous dS and dNdS values per cluster
    if ( $dS == 0 || $dS eq 'nan' ) {
        $errors{$cluster}->{'dS'}++;
        next;
    }
    if ( $dNdS == -1 || $dNdS == 99 ) {
        $errors{$cluster}->{'dNdS'}++;
        next;
    }

    # create new data structure if not already present
    unless ( $data{$cluster} ) {
        my $hash = {
            'seq1' => [],
            'seq2' => [],
            'dNdS' => [],
            'dN'   => [],
            'dS'   => [],
        };
        $data{$cluster} = $hash;
    }

    # write data
    my $i = 1;
    foreach my $key ( 'seq1', 'seq2', 'dNdS', 'dN', 'dS' ) {
        push @{ $data{$cluster}->{$key} }, $fields[ $i++ ];
    }
}

# calculate min, max, median, lower quartile, upper quartile, stdev
my @cluster_stats;   # store statistics here so we can sort them before printing
foreach my $cluster ( keys %data ) {
    my @dNdS_values = @{ $data{$cluster}->{'dNdS'} };

    # create stats object and load data
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@dNdS_values);

    # calculate stats
    my @dNdS_quantiles;
    foreach my $i ( 0 .. 4 ) {
        $dNdS_quantiles[$i] = $stat->quantile($i);
    }
    my $stdev = $stat->standard_deviation();

    my $pair_count = scalar @dNdS_values;
    my $stats = [ $cluster, $pair_count, @dNdS_quantiles, $stdev ];
    push @cluster_stats, $stats;
}

# print header
my $join_char     = $tab ? "\t" : ' ';
my @header_codes  = qw(%10s %6s %6s %6s %6s %6s %6s %6s %6s %6s);
my $header_string = join( $join_char, @header_codes ) . "\n";
printf( $header_string,
    'cluster', 'size', 'min',   '1qtr',  'med',
    '3qtr',    'max',  'stdev', 'dSerr', 'dNdSerr' );

# choose sort method
my $sort_method;
if ( $sort && $sort eq 'median' ) {
    $sort_method = \&by_median;
}
elsif ( $sort && $sort eq 'cluster' ) {
    $sort_method = \&by_cluster;
}
else { $sort_method = \&by_stdev_then_median; }

foreach my $cluster ( sort $sort_method @cluster_stats ) {

    # grab error values
    my $clust_name = $cluster->[0];
    my $dSerr =
      defined $errors{$clust_name}->{'dS'} ? $errors{$clust_name}->{'dS'} : 0;
    my $dNdSerr =
      defined $errors{$clust_name}->{'dNdS'}
      ? $errors{$clust_name}->{'dNdS'}
      : 0;

    # one code for each column
    my @printf_codes = qw(%10s %6d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6d %6d);

    # median
    if ( $cluster->[4] < 0.5 ) {
        push @printf_codes, '<';
    }

    # median >> 1 (positive selection)
    elsif ( $cluster->[4] > 1.5 ) {
        push @printf_codes, '>';
    }

    # median approx 1
    elsif ( $cluster->[4] > 0.8 && $cluster->[4] < 1.2 ) {
        push @printf_codes, '=';
    }

    # somewhere in between
    else {
    }

    my $printf_string = join( $join_char, @printf_codes ) . "\n";
    printf( $printf_string, @$cluster, $dSerr, $dNdSerr );
}

sub by_cluster {
    my ($a_number) = $a->[0] =~ /cluster(\d+)/;
    my ($b_number) = $b->[0] =~ /cluster(\d+)/;

    $a_number <=> $b_number;
}

sub by_median {
    $a->[4] <=> $b->[4];
}

sub by_stdev_then_median {
    $a->[7] <=> $b->[7]
      || $a->[4] <=> $b->[4];
}
