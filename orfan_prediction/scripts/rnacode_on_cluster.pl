#!/usr/bin/env perl
# $Id: rnacodeOnCluster.pl 2000 2010-06-09 09:34:55Z dmessina $

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Basename;

my ($dir, $outdirpath, $jobspernode, $sequential, $debug, $dryrun, $template);
GetOptions( 'outdir:s'      => \$outdirpath,
	    'indir:s'       => \$dir,
            'jobspernode:i' => \$jobspernode,
	    'sequential'    => \$sequential,
            'debug'         => \$debug,
            'dryrun'        => \$dryrun,
	    'template:s'      => \$template,
            );

my $usage = "
rnacode_on_cluster.pl - run RNAcode on a batch of alignments on the cluster

Usage: rnacode_on_cluster.pl --indir <dirpath> --outdir <dirpath> --template <file_template>

where dir_of_seqs is a directory with multiple alignment sequence files in it.
NOTE! sequence files must be in clustalw format and end in .aln

--indir          where input aligments are
--outdir         where output data should be written
--jobspernode    how many alignments to run per node
--dryrun         don't actually submit the jobs; just see the commands

--template       give a model filename with the location of the number as ?
                 e.g. myseq100.aln would be myseq?.aln

";
die $usage unless defined $dir && defined $outdirpath && $template;
die "no indirpath!\n$usage\n" unless defined $dir;
$dir = File::Spec->rel2abs($dir);
if (! -e $dir) {
    die "$dir does not exist!\n";
}
die "no outdirpath!\n$usage\n" unless defined $outdirpath;
$outdirpath = File::Spec->rel2abs($outdirpath);
if (! -e $outdirpath) {
    die "$outdirpath does not exist!\n";
}

# Globals
my $bindir = '/afs/pdc.kth.se/home/d/dmessina/bin';
my $exe   = 'RNAcode_x64';
my $exepath = File::Spec->catfile($bindir, $exe);

my $runTime = 2000;


# count files
my @files;
opendir( my $dirhandle, $dir ) or die "couldn't open $dir:$!\n";
while ( my $entry = readdir($dirhandle) ) {
    if ( $entry =~ /\.aln$/ ) {
	push @files, $entry;
    }
}
closedir($dirhandle);
my $count = @files;

# extract template
$template =~ s/\?/\$f/;

# main loop

my $node = 0;
for (my $i=0; $i < $count; $i += $jobspernode) {
    my $start = $i;
    my $end   = ($i+$jobspernode-1 >= $count) ? $count : $i+$jobspernode-1;

    my $filename = $node++ . '.sh';
    my $shscript;
    unless (defined $dryrun) {
	open($shscript, '>', $filename) or die "couldn't open $filename:!$\n";
	print $shscript '#!/bin/sh', "\n";
    }

    my ( $base, $path, $suffix ) = fileparse( $template, qr/\.[^.]*/ );
    my $outbase = $base . '.rnacode';
    my $outfile = File::Spec->catfile($outdirpath, $outbase);
    my $epsbase = $base . '_eps';
    my $epsdir  = File::Spec->catfile($outdirpath, $epsbase);
    my $errbase = $base . '.err';
    my $errfile = File::Spec->catfile($outdirpath, $errbase);
    my $infile  = File::Spec->catfile($dir, $template);
    my $args = "--tabular --outfile $outfile --eps --eps-dir $epsdir ";
    $args = $args . $infile . " " . " 2> $errfile ";

    my $command = 'for f in {'. $start . '..' . $end . '} ; do ' . $exepath . ' ' . $args . ' ; done';

    if ($dryrun) {
	print STDERR $command, "\n";
    }
    else {
	print $shscript $command, "\n";
	close($shscript);
    }

    # submit job
    my $esubmit = "esubmit -v -n 1 -t $runTime /bin/bash $filename";

    if ($dryrun) {
	print STDERR $esubmit, "\n";
    }
    else {
	my $retval = system($esubmit);

	if ( $retval > 0 ) { print STDERR $retval, "\t$!\n"; }
	else {
	    print STDERR ( "submitted job ", $node, "\n" );
	}
    }
}
