#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-06-17

use Modern::Perl;
use Getopt::Long;
use File::Basename;
use File::Spec;
use Bio::AlignIO;
my $alndir = '';
my $outdir = '';
GetOptions('alndir:s' => \$alndir,
           'outdir:s' => \$outdir,
           );

my $usage = "extract_blocks.pl - copy contiguous parts of an MSA into new MSAs

Usage: extract_blocks.pl coords_file aln_file

Options:
--alndir <dir>   location of aligned FASTA files (default: current dir)
--outdir <dir>   destination of new clustalw .aln files (default: current dir)

";
die $usage unless @ARGV == 2;

open(my $fh, '<', $ARGV[0]) or die "couldn't open $ARGV[0]:$!\n";
my $aln_file = $ARGV[1];

my $blocks;
while (my $line = <$fh>) {
    chomp $line;
    push @$blocks, $line;
}
close $fh;


my $aln_path = File::Spec->catfile($alndir, $aln_file);
my $aln_io = Bio::AlignIO->new(
    '-file'   => $aln_path,
    '-format' => 'fasta');
my $aln = $aln_io->next_aln;
unless (defined $aln) {
    die "couldn't get alignment from $aln_path\n";
}
my $aln_blocks = excise_blocks($blocks, $aln);    
my ( $base, $path, $suffix ) = fileparse( $aln_path, qr/\.[^.]*/ );
my $count = 0;
foreach my $aln_block (@$aln_blocks) {
    my $filename = $base . '_' . $count++ . '.aln';
    my $outfile = File::Spec->catfile($outdir, $filename);
    my $out = Bio::AlignIO->new(
        '-file'   => ">$outfile",
        '-format' => 'clustalw',
        '-displayname_flat' => 1,
    );
    $out->write_aln($aln_block);
}


sub process_flanks_string {
    my ($flanks_string, $aln) = @_;
    my @blocks = ($flanks_string =~ /\[(\d+\s+\d+)\]/g);
    return \@blocks;
}

sub excise_blocks {
    my ($blocks, $aln) = @_;
    
    my @aln_blocks;
    foreach my $block (@$blocks) {
        my ($start, $end) = split /\s+/, $block;
        if (defined $start && defined $end) {
            my $aln_block = $aln->slice($start, $end);
            push @aln_blocks, $aln_block;
        }
        else {
            die "malformed block coordinates: $block\n";
        }
    }
    
    return \@aln_blocks;
}