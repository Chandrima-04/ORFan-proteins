#!/usr/bin/perl
# $Id: cluster_compare_exact.pl 2947 2011-04-26 13:32:41Z dmessina $
# created by Dave Messina on 2011-04-26

use Modern::Perl;
use Cluster;
use ClusterGroup;

my $usage = "
cluster_compare_exact.pl - identify clusters which match exactly
                           between two clusterings of same data.
                           
Usage: cluster_compare_exact.pl <one.cl> <two.cl>
";
die $usage unless @ARGV == 2;

# set up clustergroup 1 as truth
my $cluster1_group  = ClusterGroup->new( '-load' => $ARGV[0] );
my @cluster1_list   = $cluster1_group->get_all_clusters(-sorted => 1);
my $cluster1_count  = scalar @cluster1_list;
my %cluster1_lookup;

foreach my $cluster (@cluster1_list) {
    my @member_ids = sort keys %{ $cluster->get_all_members };
    my $id_string = join '', @member_ids;
    $cluster1_lookup{$id_string} = $cluster->id();
}

# set up clustergroup 2
my $cluster2_group  = ClusterGroup->new( '-load' => $ARGV[1] );
my @cluster2_list   = $cluster2_group->get_all_clusters(-sorted => 1);

my $cluster2_count  = scalar @cluster2_list;
my $identical_cl_count;   # how many identical clusters?
my $identical_seq_count;  # how many sequences in those identical clusters?
my $total_seq_count; # how many sequences in total?
# check each cluster in 2 against 1
#print STDOUT join("\t", 'id', 'size'), "\n";
foreach my $cluster (@cluster2_list) {
    my @member_ids = sort keys %{ $cluster->get_all_members };
    my $id_string = join '', @member_ids;
    $total_seq_count += $cluster->size;
    
    if (exists $cluster1_lookup{$id_string}) {
        $identical_cl_count++;
        $identical_seq_count += $cluster->size;
#        print STDOUT join ("\t", $cluster->id, $cluster->size), "\n";
    }
}

my $pc_of_1    = sprintf("%d", $identical_cl_count/$cluster1_count * 100);
my $pc_of_2    = sprintf("%d", $identical_cl_count/$cluster2_count * 100);
my $pc_of_seqs = sprintf("%d", $identical_seq_count/$total_seq_count * 100);
print STDOUT "$identical_cl_count identical clusters (",
    $pc_of_1, '% of clustering 1, ', $pc_of_2, '% of clustering 2),',
    "\ncomprised of $identical_seq_count sequences (", $pc_of_seqs, "%).\n";