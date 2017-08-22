#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-05-20

use Modern::Perl;
use Getopt::Long;
use File::Basename;
use File::Spec;
use Bio::AlignIO;

my $minseqs;
my $format = 'clustalw';
GetOptions("minseqs:i" => \$minseqs,
           "format:s"  => \$format);

my $usage = "alnblocks.pl - break a multiple alignment into blocks

alnblocks.pl --minseqs 10 file1.aln file2.aln

--minseqs <n>  minimum number of seqs required for a block
--format       input alignment format (default:clustal)
";
die $usage unless @ARGV && defined $minseqs;

foreach my $file (@ARGV) {
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    
    my $in = Bio::AlignIO->new(
        '-file'   => $file,
        '-format' => $format
    );
    
    while ( my $aln = $in->next_aln() ) {
        my ($aln_by_column, $seq_lookup, $col_gaps) = pivot_aln($aln);
        find_blocks($col_gaps, $aln);
        1;        
    }
}

# given a SimpleAlign object, put sequences into an array of arrays,
# where the outer array corresponds to a column in the alignment
# and the inner array corresponds to each sequence in the alignment
# so
# seq1 AGTC
# seq2 GGGG
# seq3 CCCC
# would be
# [0] = [0] A [1] G [2] C
# [1] = [0] G [1] G [2] C
# [2] = [0] T [1] G [2] C
# [3] = [0] C [1] G [2] C
sub pivot_aln {
my ($aln) = @_;

my @aln_by_column;
my @col_gaps;
my %seq_lookup; # array index of each sequence
my $aln_length = $aln->length;

# flip alignment into array of columns
my $seq_count = 0;
foreach my $seq ($aln->each_seq) {
    $seq_lookup{$seq->id} = $seq_count;
    my $seqstring = $seq->seq;
    my @seq_array = split '', $seqstring;
    for (my $pos = 0; $pos < @seq_array; $pos++) {
        my $residue = $seq_array[$pos];
        $aln_by_column[$pos][$seq_count] = $residue;

        if ($residue eq '-') { $col_gaps[$pos]++; }
        elsif (!exists $col_gaps[$pos]) { $col_gaps[$pos] = 0; }
    }

    $seq_count++;
}

# check length of new arrays
if (scalar @aln_by_column != $aln_length) {
    die "aln_by_column has " . scalar @aln_by_column .
        "residues, but the original alignment had " . $aln_length . "\n";
}
if (scalar @col_gaps != $aln_length) {
    die "col_gaps has " . scalar @col_gaps .
        "residues, but the original alignment had " . $aln_length . "\n";
}

    return (\@aln_by_column, \%seq_lookup, \@col_gaps);
}


sub find_blocks {
    my ($col_gaps, $aln) = @_;

    my $maxgaps              = $aln->num_sequences - $minseqs;
    my %blocks;
    my $min_block_length     = 10;
    my $contiguous_gap_count = 0;
    my $current_block_length = 0;
    my %current_block;
    my $in_block
    for (my $i=0; $i < @$col_gaps; $i++) {

        if ($col_gaps->[$i] >= $maxgaps) {
            $contiguous_gap_count++;
        }
        else {
            $current_block_length++;
        }
        
        # we've hit our threshold and are starting a block
        if ($current_block_length >= $min_block_length) {
            $contiguous_gap_count = 0;
            $current_block{'start'} = $i-$current_block_length;
            print STDERR 'block started at ' . $current_block{'start'} . "\n";
        }
        
    }
}

#    my $i = 0;
# my $outbase = $base . '_' . $i++ . '.aln';
# my $outfile = File::Spec->catfile( $path, $outbase );

