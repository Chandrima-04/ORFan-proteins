#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-04-28

use Modern::Perl;
use Getopt::Long;
use Bio::AlignIO;

# my ();
# GetOptions();

my $usage = "
alns2dist.pl - reduce multiple alignments to have the number of seqs
               specified by a distfile

usage: alns2dist.pl <stockholm aln file> <distfile>

where
distfile is a space-delimited file with two columns: # of seqs, and count
e.g. 1 255
means 255 alignments of 1 sequence each.
";
die $usage unless @ARGV;

my %dist;

open my $dist, '<', $ARGV[1] or die "couldn't open $ARGV[1]\n";
while (my $line = <$dist>) {
    chomp $line;
    my ($num, $count) = split /\s+/, $line;
    next if $count == 0;
    $dist{$num} = $count;
}

my $aln_in = Bio::AlignIO->new(-format => 'stockholm',
                               -file   => $ARGV[0] );

# my $aln_out = Bio::AlignIO->new(-format => 'fasta',
#                                 -fh     => \*STDOUT );

while (my $aln = $aln_in->next_aln) {
    my $seq_count   = $aln->num_sequences;
    my $num_to_keep = num_needed($seq_count);

    # we've filled our distribution and are all done
    if (!defined $num_to_keep) { last; }
    else {
        my $new_aln     = $aln->select(1, $num_to_keep);
        foreach my $seq ($new_aln->each_seq) {
            my $id = $seq->id;
            if (defined $seq->{'_version'}) {
                $id = $id . '.' . $seq->{'_version'};
            }
            print STDOUT join("\t", $aln->id, $id), "\n";
        }
    }
}

sub num_needed {
    my ($seq_count) = @_;
    my $needed_count;

    # start with the largest # of seqs needed and go down
    my @sorted_dist = sort { $b <=> $a } keys %dist;
    foreach my $count (@sorted_dist) {
        if ($seq_count >= $count) {
            $needed_count = $count;
            if ($dist{$count} == 1) {
                delete $dist{$count};
            }
            else { $dist{$count}--; }
            
            return $needed_count;
        }
    }
    
    if (!defined $needed_count) {
        return;
    }
}