#!/usr/bin/perl
# $Id: batch_compare_clusters.pl 2935 2011-04-19 09:49:49Z dmessina $
# created by Dave Messina on 2011-04-19

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($n);
GetOptions("n=i" => \$n);

my $usage = "
batch_compare_clusters.pl --n <number of jobs> one.cl two.cl
";
die $usage unless (@ARGV == 2 && defined $n);

# grab command-line args
my $one = $ARGV[0];
my $two = $ARGV[1];
my ( $base, $path, $suffix ) = fileparse( $two, qr/\.[^.]*/ );

# get number of clusters in the 2nd cluster file
my $stat_raw = `cluster_stat $two`;
my $total_clusters;
if ($stat_raw =~ /(\d+)\sclusters/) {
    $total_clusters = $1;
}

# determine clusters per jobs
my $clusters_per_job = sprintf("%d", $total_clusters/$n);
my $remainder = $total_clusters % $n;
if ( ($clusters_per_job * $n) + $remainder != $total_clusters ) {
    die "won't process all the jobs!\n";
}

# submit jobs
print STDERR 'total: ', $total_clusters, "\n";
for (my $i = 0; $i < $total_clusters; $i += $clusters_per_job) {
    my $start = $i;
    my $end   = $i+$clusters_per_job-1;
    if ($end >= $total_clusters) {
        $end = $total_clusters;
    }
    my $out = $path . $i . '.bccout';
    print STDERR join("\t", $start, $end, $out), "\n";
    `esubmit -n 1 -t 2000 compare_clusters.pl --start $start --end $end --print --out $out $one $two`
}
