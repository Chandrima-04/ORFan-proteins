#!/usr/bin/perl
# $Id: cluster_string2stats.pl 2940 2011-04-19 14:16:48Z dmessina $
# created by Dave Messina on 2011-04-19

use strict;
use warnings;
use ClusterResult;

my $usage = "cluster_string2stats.pl <strings as files>

compiles output from separate runs of compare_clusters.pl
into summary stats.
";

my ($ss, $sd, $ds, $dd);
my $result = ClusterResult->new();

while (my $line = <>) {
    chomp $line;
    
    my @fields = split "\t", $line;
    if (scalar @fields == 6) {
        $ss += $fields[2];
        $sd += $fields[3];
        $dd += $fields[4];
        $ds += $fields[5];
    }
    elsif ($fields[0] eq 'SD') { 
        $result->add_seq_to_SD('seqid'   => $fields[1],
                               'cluster' => $fields[2]);
    }
    elsif ($fields[0] eq 'DS') {
        $result->add_seq_to_DS('seqid'   => $fields[1],
                               'cluster' => $fields[2]);        
    }
    else {
        die "unexpected line $line\n";
    }
}

$result->SS($ss);
$result->SD($sd);
$result->DD($dd);
$result->DS($ds);
unless ( $result->is_corrected ) {
    $result->correct;
}

$result->tabulate_SD_and_DS;


print 'SS:         ', $result->SS, "\n";
print 'DD:         ', $result->DD, "\n";
print 'SD          ', $result->SD, "\n";
print 'DS:         ', $result->DS, "\n";
print 'rand:       ', $result->rand, "\n";
print 'jaccard:    ', $result->jaccard(), "\n";
print 'fmi:        ', $result->folkes_mallows(), "\n";
print 'in cluster: ', $result->in_cluster(), "\n";
print 'seq diff:   ', $result->seq_diff(), "\n";