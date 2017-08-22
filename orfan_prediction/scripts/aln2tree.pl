#!/usr/bin/perl
# $Id: aln2tree.pl 1607 2010-02-12 14:26:15Z dmessina $

use strict;
use warnings;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;
use Bio::Tree::DistanceFactory;
use Bio::Matrix::IO;
use Bio::TreeIO;
use File::Basename;
use Getopt::Long;

my ($dist);
GetOptions('dist:s' => \$dist);

my $usage =
"aln2tree.pl - takes multiple alignments in FASTA and makes NJ trees
               using the Kimura distance method.
               
Usage: aln2tree.pl aln1 .fa [aln2.fa] ... [aln3.fa]

Output filename(s) are based on input filenames, changing the suffix to .tree.

Options:
--dist  <suffix>   take distances from external file(s) in phylip distance
                   format. Files should be named same as input files but
                   with the <suffix> attached.
                   
                   e.g. if --dist jc
                   cluster132_dnatrim.fasta
                   cluster132_dnatrim.fasta.jc
";
die $usage unless @ARGV;

while ( my $file = shift @ARGV ) {
    my ( $name, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my $outname = $path . '/' . $name . '.tree';
    
    my $alnio = Bio::AlignIO->new( -file => $file, -format => 'fasta' );
    my $dfactory = Bio::Tree::DistanceFactory->new( -method => 'NJ' );
    my $stats    = Bio::Align::DNAStatistics->new;
    my $treeout  = Bio::TreeIO->new( -format => 'newick',
                                     -file   => ">$outname");
    while ( my $aln = $alnio->next_aln ) {
        my $matrix;
        if ($dist) {
            my $distfile = $file . '.' . $dist;
            die ("no distfile $distfile") if (!-e $distfile);
            my $parser = Bio::Matrix::IO->new(-format => 'phylip',
                                              -file   => $distfile);
            $matrix = $parser->next_matrix;
        }
        else {
            $matrix = $stats->D_Kimura_variance($aln);
        }
        my $tree = $dfactory->make_tree($matrix);
        $treeout->write_tree($tree);
    }
}
