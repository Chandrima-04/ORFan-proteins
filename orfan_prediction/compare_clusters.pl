#!/usr/bin/perl
# $Id: compare_clusters.pl 2969 2011-05-11 11:10:53Z dmessina $
# created by Dave Messina on 2011-04-14

use strict;
use warnings;
use Cluster;
use ClusterGroup;
use ClusterUtils;
use ClusterResult;
use Getopt::Long;
# use Data::Dumper;
# $Data::Dumper::Sortkeys = 1;

my ($cl_start, $cl_end, $print, $verbose, $outfile);
GetOptions(
    "start=i" => \$cl_start,
    "end=i"   => \$cl_end,
    "print"   => \$print,
    "verbose" => \$verbose,
    "out=s"   => \$outfile,
);

die usage() unless @ARGV == 2;

my $outfh;
if (defined $outfile) {
    open ($outfh, '>', $outfile) or die "couldn't open $outfile: $!\n";
}
else {
    $outfh = \*STDOUT;
}


# initialize result object
my $result = ClusterResult->new();

# set up clustergroup 1 as truth
my $cluster1_group  = ClusterGroup->new( '-load' => $ARGV[0] );
my @cluster1_list   = $cluster1_group->get_all_clusters;
my $cluster1_matrix = {};
my @cluster1_allseqids;

# populate true matrix
foreach my $cluster (@cluster1_list) {
    my @member_ids = keys %{ $cluster->get_all_members };
    push @cluster1_allseqids, @member_ids;
    add_cluster_to_matrix( \@member_ids, $cluster1_matrix );
}
print STDERR "built true matrix\n" if $verbose;

#print Dumper $cluster1_matrix;

# set up clustergroup 2
my $cluster2_group  = ClusterGroup->new( '-load' => $ARGV[1] );
my @cluster2_list   = $cluster2_group->get_all_clusters;
my $cluster2_matrix = {};

# check and process start and end
$cl_start = defined $cl_start ? $cl_start : 0;
$cl_end   = defined $cl_end   ? $cl_end   : $#cluster2_list;
if ($cl_start < 0 || $cl_end > $#cluster2_list) {
    die "start or end is out of bounds!"
}
$result->cluster_start($cl_start);
$result->cluster_end($cl_end);

# counter for user feedback
# my $total_clusters  = scalar @cluster2_list;
# my $one_percent     = $total_clusters * .01;
# if ($one_percent < 1) { $one_percent = 1; }
# my $current_cluster = 0;

# test clustergroup 2 against clustergroup 1
for (my $i = $cl_start; $i <= $cl_end; $i++) {
    my $cluster2 = $cluster2_list[$i];
    test_cluster(
        $cluster2,        \@cluster1_allseqids,
        $cluster2_matrix, $cluster1_matrix
    );

    # $current_cluster++;
    # if ( $verbose && $current_cluster % $one_percent == 0 ) {
    #     printf STDERR "%d%%\n", ($current_cluster/$total_clusters) * 100;
    # }
}

if ($print) {
    print $outfh $result->to_string, "\n";
    while (my ($key, $value) = each %{$result->SD_seqs}) {
        print $outfh join("\t", 'SD', $key, $value), "\n";
    }
    while (my ($key, $value) = each %{$result->DS_seqs}) {
        print $outfh join("\t", 'DS', $key, $value), "\n";
    }
}

if ($verbose) {
    print "\n";
    print 'rand:       ', $result->rand, "\n";
    print 'jaccard:    ', $result->jaccard(), "\n";
    print 'fmi:        ', $result->folkes_mallows(), "\n";
    print 'in cluster: ', $result->in_cluster(), "\n";
    print 'seq diff:   ', $result->seq_diff(), "\n";
}

#print STDERR Dumper $cluster2_matrix;
# print Dumper \@cluster1_allseqids;

sub usage {
    print STDERR "
compare_clusters.pl [options] one.cl two.cl

Options:
--start             id of cluster to start with (default 0)
--end               id of cluster to end with   (default last cluster)
--print             print a one-line output of result numbers
--out <filename>    name of file to write output to (default STDOUT)
--verbose           give more feedback

";
    
    exit 1;
}

sub add_cluster_to_matrix {
    my ( $seq_ids, $matrix ) = @_;

    for ( my $i = 0 ; $i < @$seq_ids ; $i++ ) {
        for ( my $j = 0 ; $j < @$seq_ids ; $j++ ) {
            my $seq_i = $seq_ids->[$i];
            my $seq_j = $seq_ids->[$j];
            next if $seq_i eq $seq_j;
            $matrix->{$seq_i}->{$seq_j} = 'SS';
        }
    }
}

# populate test_matrix with the results of comparing the seq ids from a
# cluster to another cluster group's matrix (true_matrix)
sub test_cluster {
    my ( $cluster, $seq_ids, $test_matrix, $true_matrix ) = @_;

    my @member_ids = keys %{ $cluster->get_all_members };

    # first check each seqid against the members of the same cluster
    for ( my $i = 0 ; $i < @member_ids ; $i++ ) {
        for ( my $j = 0 ; $j < @member_ids ; $j++ ) {
            my $seq_i = $member_ids[$i];
            my $seq_j = $member_ids[$j];
            next if $seq_i eq $seq_j;

            my $value_in_true = $true_matrix->{$seq_i}->{$seq_j};
            if ( defined $value_in_true && $value_in_true eq 'SS' ) {
                $test_matrix->{$seq_i}->{$seq_j} = 'SS';
                $result->SS(1);
            }
            else {
                $test_matrix->{$seq_i}->{$seq_j} = 'DS';
                $result->DS(1);
                $result->add_seq_to_DS('seqid'   => $seq_j,
                                       'cluster' => $cluster->id);
            }
        }
    }

    # make a lookup of this cluster's members so they aren't checked again
    my %this_cluster = map { $_ => 1 } @member_ids;

    # then check each seqid against all of the other remaining seqids
    for ( my $i = 0 ; $i < @member_ids ; $i++ ) {
        for ( my $j = 0 ; $j < @$seq_ids ; $j++ ) {
            my $seq_i = $member_ids[$i];
            my $seq_j = $seq_ids->[$j];
            next if $this_cluster{$seq_j};
            next if $seq_i eq $seq_j;

            my $value_in_true = $true_matrix->{$seq_i}->{$seq_j};
            if ( defined $value_in_true && $value_in_true eq 'SS' ) {
                $test_matrix->{$seq_i}->{$seq_j} = 'SD';
                $result->SD(1);
                $result->add_seq_to_SD('seqid'   => $seq_j,
                                       'cluster' => $cluster->id);
            }
            else {
                $test_matrix->{$seq_i}->{$seq_j} = 'DD';
                $result->DD(1);
            }
        }
    }
}
