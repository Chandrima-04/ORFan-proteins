#!/usr/bin/perl
# $Id: remap_codeml_data.pl 1666 2010-03-15 17:09:39Z dmessina $

use strict;
use warnings;
use File::Find;
use Env qw(PWD);

my $usage = "
remap_codeml_data.pl cluster_seq_pairs pairs_codeml new_output_dir

where cluster_seq_pairs is a directory containing subdirectories, one for
each cluster:

cluster_seq_pairs/cluster220_dnatrim_uniq/
                 /cluster221_dnatrim_uniq/
                 ...

and in each of those cluster dirs are the pair fasta files

cluster_seq_pairs/cluster220_dnatrim_uniq/all_c1088_all_c48754.fa
                                         /all_c31507_all_c1088.fa
                                          ...
                                          
and where pairs_codeml is a directory containing subdirectories, one for
each cluster:

pairs_codeml/cluster220/
            /cluster221/
            ...

and in each of those cluster dirs are randomly-named subdirectories which
contain the codeml output files:

pairs_codeml/cluster220/davepaml614/all_c1088_all_c48754.out
                                   /all_c31507_all_c1088.out
                                   ...

and where new_output_dir is the name of an (already existing) directory
where the new cluster dirs will be created.

This program is designed to cleanup an error by which some codeml output files
are put in the wrong cluster directory.

It will make a lookup table of all the pairs from the
cluster_seq_pairs fasta files with their cluster name, and use it to copy
the pairs_codeml output files to a new set of cluster directories.

It will also print out an error file called 'remap.errors', which has 3
columns: the pair, its real cluster, and its incorrect cluster.
";
die $usage unless @ARGV == 3;

my ( $seqdir, $codemldir, $outdir ) = @ARGV;
my %pair_cluster;    # maps seq pairs to a cluster name

open( ERR, ">remap.errors" ) or die "coudln't create remap.errors file:$!\n";

# build lookup hash
find( \&make_fasta_cluster_hash, $seqdir );

# mv codeml output files into new directoris
find( \&move_codeml_output, $codemldir );

sub make_fasta_cluster_hash {
    return unless -f;               # only look at files
    return unless $_ =~ /\.fa$/;    # filename must end in .fa

    # get cluster name
    my $cluster = get_cluster_name($File::Find::dir);

    # get pair name
    my $pairname = $_;
    $pairname =~ s/\.fa//;

    if ( $pair_cluster{$pairname} ) {
        print STDERR "WARNING: duplicate:$pairname was:",
          $pair_cluster{$pairname}, " is:", $cluster, "\n";
    }
    else {
        $pair_cluster{$pairname} = $cluster;
    }
}

sub move_codeml_output {
    return unless -f;                # only look at files
    return unless $_ =~ /\.out$/;    # filename must end in .out

    # extract cluster name
    my $cluster = get_cluster_name($File::Find::dir);

    # get pair name
    my $filename = $_;
    my $pairname = $_;
    $pairname =~ s/\.out//;

    if ( $pair_cluster{$pairname} ) {
        my $real_cluster = $pair_cluster{$pairname};

        # print cluster mismatches to errfile
        if ( $cluster ne $real_cluster ) {
            print ERR join( "\t", $pairname, $real_cluster, $cluster ), "\n";
        }

        # make output cluster dir unless it exists
        my $outpath = join('/', $ENV{'PWD'}, $outdir, $real_cluster);
        unless ( -d $outpath ) {
            mkdir($outpath, 0777) or die "couldn't make dir $outpath:$!\n";
        }

        # copy codeml outfile to outdir
        system("cp $filename $outpath")
          && die "couldn't copy $filename to $outpath:$!\n";
    }
    else {
        die "couldn't find $pairname in lookup!\n";
    }

}

sub get_cluster_name {
    my ($filename) = @_;

    # extract cluster name
    my $cluster;
    if ( $filename =~ /(cluster\d+)/ ) {
        $cluster = $1;
    }
    else {
        die "couldn't find clustername in $filename\n";
    }

    return $cluster;
}
