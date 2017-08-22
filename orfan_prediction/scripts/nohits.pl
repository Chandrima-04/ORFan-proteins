#!/usr/bin/perl
# $Id: nohits.pl 1819 2010-05-05 11:37:49Z dmessina $
# created by Dave Messina on 2010-05-05

use Modern::Perl;
use File::Basename;

my $usage = "nohits.pl - identify BLAST queries that had no hits
    Usage: nohits.pl <NCBI-BLAST tabular output files>

    output is the query IDs of the no-hitters.
";
die $usage unless @ARGV;


while ( my $file = shift @ARGV ) {
    open( my $fh, $file ) or die "couldn't open $file: $!\n";
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    die "no filename base found in $file!\n" unless defined $base;

    my $querycount   = 0;
    my $hitlinecount = 0;
    while ( my $line = <$fh> ) {
        chomp $line;

        # new record
        if ( $line eq '# BLASTX 2.2.23+' ) {
            my ( $query, $db, $hitcount );

            # query
            my $queryline = <$fh>;
            chomp $queryline;
            if ( $queryline =~ /^# Query:\s+(.+?)\s\/accession/ ) {
                $query = $1;
                $querycount++;
            }
            else { die "couldn't read queryline on $base line $.: $queryline\n"; }

            # db
            my $dbline = <$fh>;
            chomp $dbline;
            if ($dbline =~ /^# RID:/) {
                $dbline = <$fh>;
                chomp $dbline;
            }
            if ( $dbline =~ /^# Database:\s+(.+)$/ ) {
                $db = $1;
            }
            else { die "couldn't read dbline on $base line $.: $dbline\n"; }

            # hits
            my $hitline = <$fh>;
            chomp $hitline;
            if ( $hitline =~ /^# Fields:/ ) {
                $hitline = <$fh>;
                chomp $hitline;
            }
            if ( $hitline =~ /^# (\d+) hits found$/ ) {
                $hitcount = $1;
                $hitlinecount++;
                if ( $hitcount == 0 ) {
                    $hitcount = '-none-';
                    say $query;

                    #say join("\t", $query, $db, $hitcount);
                }
            }
            else { die "couldn't read hitline: $hitline  on $base line $.\n"; }
        }
        elsif ( $line =~ /^# BLAST processed (\d+) quer/ ) {
            my $official_count = $1;
            if (   $querycount != $official_count
                || $hitlinecount != $official_count )
            {
                die "count discrepancy in $base! query: $querycount ",
                  "hit: $hitlinecount official: $official_count\n";
            }
            else {
                $querycount   = 0;
                $hitlinecount = 0;
            }
        }
        else { next; }
    }
}