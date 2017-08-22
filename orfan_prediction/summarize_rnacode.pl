#!/usr/bin/perl
# $Id: summarize_rnacode.pl 2695 2010-12-07 10:14:16Z dmessina $
# created by Dave Messina on 2010-05-28

use Modern::Perl;
use File::Basename;
use Statistics::Descriptive;

my $usage = "
summarize_rnacode.pl - summarize the results from multiple RNAcode runs

Usage: summarize_rnacode.pl <lots of rnacode output files>

Note that the rnacode output needs to be in tab format and named by clusters

e.g. cluster139_fake_10.rnacode

";
die $usage unless @ARGV;

# global storage
my %cluster_best_pval; # hash of arrays

FILE: foreach my $file (@ARGV) {

    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    
    # extract cluster name from path
    # need to take the regexp from a command-line parameter
    my ($cluster, $replicate) = $base =~ /(cluster\d+)_fake_(\d+)/;
#    my ($cluster, $replicate) = $base =~ /(methjann-(?:non)*coding_\d+)_(\d+)/;
    die "no cluster name in $path\n" unless $cluster;
    die "no replicate number in $path\n" unless $replicate;

    # check for empty file and give it a high pval
    # if (-z $file) {
    #     push @{ $cluster_best_pval{$cluster} }, 1;
    #     next FILE;
    # }

    open (my $fh, '<', $file) or die "couldn't open $file: $!\n";
    LINE: while (my $line = <$fh>) {

        chomp $line;
        my @fields = split /\t/, $line;
        if (scalar @fields != 11) {
            warn "Line $. from $base doesn't look right:\n", $line, "\n";
        }
        else {
            push @{ $cluster_best_pval{$cluster} }, $fields[10];
        }
    
        last LINE; # for now we want best hit only
    }
}

# header
#printf STDOUT ( "%10s %6s %5s  %5s\n", 'Cluster', 'Median', 'Count', 'Min' );
printf STDOUT ( "%10s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\n", 'Cluster', 'Count',
    'Min', '1Q', 'Median', '3Q', 'Max' );

foreach my $cluster ( sort by_cluster keys %cluster_best_pval ) {
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data( @{ $cluster_best_pval{$cluster} } );
    # printf STDOUT ( "%10s %6.3f %5d  %6.1g\n", $cluster, $stat->median(),
    #     $stat->count(), $stat->min(), );
        
    printf STDOUT ( "%10s\t%5d\t%6.2g\t%6.2g\t%6.2g\t%6.2g\t%6.2g\n",
        $cluster, $stat->count(),
        $stat->quantile(0), $stat->quantile(1), $stat->quantile(2),
        $stat->quantile(3), $stat->quantile(4),
        );
}

sub by_cluster {
    my ($a_number) = $a =~ /cluster(\d+)/;
    my ($b_number) = $b =~ /cluster(\d+)/;
#    my ($a_number) = $a =~ /methjann-(?:non)*coding_(\d+)/;
#    my ($b_number) = $b =~ /methjann-(?:non)*coding_(\d+)/;
    $a_number <=> $b_number;
}
