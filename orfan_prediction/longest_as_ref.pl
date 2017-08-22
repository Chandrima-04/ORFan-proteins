#!/usr/bin/perl
# $Id: longest_as_ref.pl 2052 2010-06-15 11:35:14Z dmessina $
# created by Dave Messina on 2010-06-11

use Modern::Perl;
use Getopt::Long;
use File::Basename;
use File::Spec;
use Bio::AlignIO;
use File::Path qw(make_path);

my ( $outdirpath, $alnlen, $seqs, $outfmt );
GetOptions(
    'outdir:s' => \$outdirpath,
    'alnlen:i' => \$alnlen,
    'seqs:i'   => \$seqs,
    'outfmt:s' => \$outfmt,
);

my $usage = "
longest_as_ref.pl - move longest seq in multiple alignment to be first,
                    where longest here means longest ungapped sequence.

Usage: longest_as_ref.pl --outdir <outdirpath> <aln1> ... [alnN]

where
--outdir       is the directory where you want the output files to be written.
--outfmt       is the desired output format (default clusltalw)
--alnlen <n>   minimum ungapped length of an alignment to be used (skip shorties)
--seqs   <n>   minimum number of seqs in an alignment (skip small alignments)
";
die $usage unless (@ARGV);
die "no outdirpath!\n$usage\n" unless defined $outdirpath;
$outdirpath = File::Spec->rel2abs($outdirpath);
if ( !-e $outdirpath ) {
    say STDERR "creating output directory $outdirpath (didn't exist)";
    make_path($outdirpath)
      or die "couldn't create output dir $outdirpath!: $!\n";
}

foreach my $file (@ARGV) {

    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );

    my $in = Bio::AlignIO->new(
        '-file'   => $file,
        '-format' => 'clustalw'
    );

    $outfmt //= 'clustalw';
    my $outsuffix = ($outfmt eq 'clustalw') ? '.aln' : ".$outfmt";
    my $outbase = $base . '_long1st' . $outsuffix;
    my $outfile = File::Spec->catfile( $outdirpath, $outbase );
    my $out     = Bio::AlignIO->new(
        '-file'             => ">$outfile",
        '-format'           => $outfmt,
        '-displayname_flat' => 1,
    );

    while ( my $aln = $in->next_aln() ) {

        # filter if requested
        if ( defined $alnlen && $aln->num_residues < $alnlen ) {
            say STDERR "skipping aln in $file. aln too short ($alnlen res)";
            next;
        }
        if ( defined $seqs   && $aln->num_sequences < $seqs ) {
            say STDERR "skipping aln in $file. only $seqs seqs";
            next;            
        }

        # otherwise, sort and write out
        sort_by_length($aln);
        $out->write_aln($aln);
    }

}

# longest first
sub sort_by_length {
    my $aln_obj = shift;
    my ( @array, $count );
    foreach my $seq ( $aln_obj->each_seq() ) {
        push @array, $seq;
    }

    $count = 0;
    %{ $aln_obj->{'_order'} } = ();    # reset the hash;

    foreach my $seq ( sort by_length @array ) {
        my $nse = $seq->get_nse;
        my $len = $seq->end() - $seq->start() + 1;
        $aln_obj->{'_order'}->{$count} = $nse;
        print STDERR "$count\t$len\t$nse\n";
        $count++;
    }
    1;
}


sub by_length {
    my $a_len = $a->end() - $a->start() + 1;
    my $b_len = $b->end() - $b->start() + 1;
    
    $b_len <=> $a_len;
}
