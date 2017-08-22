#!/usr/bin/env perl
# $Id: batch_evolver.pl 1766 2010-04-26 13:15:17Z dmessina $

use strict;
use warnings;
use File::Basename;
use File::Spec;
use File::Copy;
use Getopt::Long;
use Bio::Tools::Run::Phylo::PAML::Evolver;
use Bio::TreeIO;

# hard-coded global - may want to change in the future
my $TREEFILEPREFIX = 'RAxML_bestTree';

my ( $codemldir, $treedir, $outdir, $numcodons );
GetOptions(
    'codemldir:s' => \$codemldir,
    'treedir:s'   => \$treedir,
    'outdir:s'    => \$outdir,
    'numcodons:i' => \$numcodons,
);

my $usage = "
batch_evolver.pl - run evolver on a bunch of data in one go

Usage: batch_evolver.pl --codemldir dir --treedir dir --outdir dir

where

codemldir is a directory containing codeml output files. These must have the
          .mlc suffix and have been run with sufficient verbosity to contain
          the codon frequency table
          
treedir   is a directory containing tree files. These files must have names
          like 'RAxML_bestTree.cluster100'.

outdir    is a directory where the evolver output files should be written.

Options:
--numcodons <num>   specify the number of codons that the generated seqs should
                    have (default: use number of codons that input seqs have)
";
die "\nno codemldir!\n$usage" unless defined $codemldir;
die "\nno treedir!\n$usage"   unless defined $treedir;
die "\nno outdir!\n$usage"    unless defined $outdir;

# check and open dirs
foreach my $dir ( $codemldir, $treedir, $outdir ) {
    die "$dir is not a directory!\n" unless -d $dir;
}
$codemldir = File::Spec->rel2abs($codemldir);
$treedir   = File::Spec->rel2abs($treedir);
$outdir    = File::Spec->rel2abs($outdir);

opendir( my $codemldirhandle, $codemldir )
  or die "couldn't open $codemldir: $!\n";

FILE: while ( my $codeml = readdir($codemldirhandle) ) {

    next unless $codeml =~ /\.mlc$/;

    # make sure we can reach that file from anywhere on the filesystem
    my $file = File::Spec->catfile( $codemldir, $codeml );
    $file = File::Spec->rel2abs($file);

    # grab cluster ID
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my $cluster;
    if ( $codeml =~ /(cluster\d+)/ ) { $cluster = $1; }
    else { die "couldn't identify cluster from filename $codeml\n"; }

    my $evolver =
      Bio::Tools::Run::Phylo::PAML::Evolver->new( '-save_tempfiles' => 1 );

    # read in codeml file
    my $parser      = Bio::Tools::Phylo::PAML->new( '-file' => $file );
    my $result      = $parser->next_result();
    if ( !defined($result) ) {
        warn "couldn't parse $file\n";
        next FILE;
    }
    my @codon_freqs = $result->get_CodonFreqs();
    $evolver->set_CodonFreqs( \@codon_freqs );

    # grab tree
    my $treefile = File::Spec->catfile( $treedir, "$TREEFILEPREFIX.$cluster" );
    die "no such file $treefile!\n" unless ( -s $treefile );
    my $treeio = Bio::TreeIO->new(
        '-file'   => $treefile,
        '-format' => 'newick'
    );
    my $tree = $treeio->next_tree();
    $evolver->tree($tree);

    # set parameters
    my @seqs = $result->get_seqs;
    my $num_of_codons = $numcodons // ($seqs[0]->length / 3);
    $evolver->set_parameter('nuclsites', $num_of_codons);
    $evolver->set_parameter( 'replicates', 100 );
    $evolver->set_parameter('tree_length', -1);
    $evolver->set_parameter('kappa', 1);
    $evolver->set_parameter('omega', 1);
    $evolver->prepare();

    # run it
    my $retval = $evolver->run();
    if ( $retval == 0 ) {    # there was a problem
        die $evolver->error_string, "\n",
          "tempdir: ", $evolver->tempdir(), "\n";
    }
    else {                   # everything's okay, grab the output file
        my $tempdir      = $evolver->tempdir();
        my $tempfile     = File::Spec->catfile( $tempdir, 'mc.paml' );
        my $outfile_name = $cluster . '_fake.phy';
        my $outfile_path = File::Spec->catfile( $outdir, $outfile_name );
        copy( $tempfile, $outfile_path )
          or die "Couldn't copy $tempfile to $outfile_path: $!\n";

        # cleanup
        $evolver->cleanup();
    }
}
