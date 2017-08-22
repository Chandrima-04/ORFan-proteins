#!/usr/bin/perl
# $Id: remove_identicals.pl 1637 2010-02-19 16:49:24Z dmessina $

use strict;
use warnings;
use Bio::AlignIO;
use Bio::SimpleAlign;
use File::Basename;

my $usage = "
remove_identicals - compare all seqs in a multi-FASTA file to each other and
                    remove any seqs which are an exact duplicate of another.
                    
                    A tab-delimited list of the unique IDs and the corresponding
                    duplicate IDs is written to a file. The filename has the
                    same basename as the corresponding input file.
                    e.g my.fa => my.dupes.txt

Usage: remove_identicals my.fa my2.fa ... myN.fa
";
die $usage unless @ARGV;

while ( my $file = shift @ARGV ) {
    my $instream = Bio::AlignIO->new(
        '-format' => 'fasta',
        '-file'   => $file,
    );

    # output sequence file
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my $outfile = $path . '/' . $base . '_uniq' . $suffix;
    my $outstream = Bio::AlignIO->new(
        '-format' => 'fasta',
        '-file'   => ">$outfile",
    );

    # output dupes file
    my $dupes_name = $path . '/' . $base . '_dupes' . '.txt';
    open(DUPES, ">$dupes_name") or die "couldn't open $dupes_name for writing. $!\n";
    
    while (my $old_aln = $instream->next_aln()) {
        my ($uniq_aln, $ids) = $old_aln->uniq_seq;
        
        my $new_aln = fix_names($uniq_aln, $ids);
        $outstream->write_aln($new_aln);
    }
    
    close(DUPES);
}


sub fix_names {
    my ($uniq_aln, $ids) = @_;
    
    my $new_aln = Bio::SimpleAlign->new();
    
    foreach my $seq ( $uniq_aln->each_seq() ) {
        my $ST_name =  $seq->display_id;

        # grab all the IDs
        my ($ST_num) = $ST_name =~ /ST(\d+)/;
        my @IDs = @{ $ids->{$ST_num} };
        die "no IDs for $ST_name!\n" unless @IDs > 0;

        # set the sequence ID to be the first listed ID
        $seq->display_id($IDs[0]);

        # add that seq to the new alignment object
        $new_aln->add_seq($seq);
        
        # write out the dupes
        print DUPES join("\t", @IDs), "\n";
    }
    
    return $new_aln;
}