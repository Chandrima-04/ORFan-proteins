#!/usr/bin/perl
# $Id: parse_multi_codeml.pl 1831 2010-05-06 12:47:49Z dmessina $

use strict;
use warnings;
use File::Basename;
use Bio::Tools::Phylo::PAML;
use IO::Scalar;

my $usage = "
parse_multi_codeml - grabs Ka/Ks info from codeml result file(s)
                     that have been run with ndata and therefore
                     have multiple results in one file and reports
                     in a tab-delimited format.

Report format is:
Cluster    Seq1    Seq2    dN/dS    dN    dS

NOTE: This program assumes the codeml result files are named after
their cluster, like 'cluster123'.

";
die $usage unless @ARGV;

while ( my $infile = shift @ARGV ) {
    my ( $name, $path, $suffix ) = fileparse( $infile, qr/\.[^.]*/ );

    # extract cluster name from path
    my ($cluster) = $name =~ /(cluster\d+)/;
    die "no cluster name in $path\n" unless $cluster;

    # slurp file, splitting it into pieces, one run per piece
    my @runs;
    {
        local $/ = undef;
        open( my $fh, $infile ) or die "couldn't open $infile";
        my $whole_file = <$fh>;

        # split if we've got multiple runs
        if ($whole_file =~ /Data set/) {
            @runs = split /\n(?=Data set \d+)/, $whole_file;

            # get rid of first piece with just the 'seed used' in it
            shift @runs;
        }
        # otherwise, assume just one run
        else {
            push @runs, $whole_file;
        }
    }

    # treat each run as a separate file
    for ( my $j = 0 ; $j < @runs ; $j++ ) {

        # cluster run ID
        # if there are multiple runs of fake seqs for each cluster
        # we need a way to identify which run and that the data is fake
        my $cluster_runID;
        if (scalar @runs > 1) {
            # grab run number to make cluster_runID
            my @lines = split "\n", $runs[$j];
          LINE: foreach my $line (@lines) {
                if ( $line =~ /^Data set (\d+)$/ ) {
                    my $run_num = $1;
                    $cluster_runID = $cluster . 'fake.' . $run_num;
                    undef @lines;
                    last LINE;
                }
            }
        }
        # otherwise, it's just a single run and we just want the cluster name
        else { $cluster_runID = $cluster; }

        # parse run
        my $run_fh   = IO::Scalar->new( \$runs[$j] );
        my $parser   = Bio::Tools::Phylo::PAML->new( '-fh' => $run_fh );
        my $result   = $parser->next_result();
        
        if (defined $result) {
            my $MLmatrix = $result->get_MLmatrix();
            my @otus     = $result->get_seqs();

            for ( my $i = 0 ; $i < ( scalar @otus - 1 ) ; $i++ ) {
                for ( my $j = $i + 1 ; $j < ( scalar @otus ) ; $j++ ) {

                    # check we've got everything
                    if (   defined($cluster_runID)
                        && defined( $otus[$i]->display_id )
                        && defined( $otus[$j]->display_id )
                        && defined( $MLmatrix->[$i]->[$j]->{'omega'} )
                        && defined( $MLmatrix->[$i]->[$j]->{'dN'} )
                        && defined( $MLmatrix->[$i]->[$j]->{'dS'} ) )
                    {
                        print STDOUT join( "\t",
                            $cluster_runID,
                            $otus[$i]->display_id,
                            $otus[$j]->display_id,
                            $MLmatrix->[$i]->[$j]->{'omega'},
                            $MLmatrix->[$i]->[$j]->{'dN'},
                            $MLmatrix->[$i]->[$j]->{'dS'},
                          ),
                          "\n";
                    }
                    else {
                        print STDERR join( "\t",
                            'ERROR:',
                            $cluster_runID,
                            $otus[$i]->display_id,
                            $otus[$j]->display_id,
                            ),
                            "\n";
                    }

                }
            }
        }
        else {
            print STDOUT join("\t", $cluster_runID, '-no sites-'), "\n";
        }
    }
}
