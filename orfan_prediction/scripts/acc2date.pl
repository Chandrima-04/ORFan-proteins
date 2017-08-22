#!/usr/bin/perl
# $Id: acc2date.pl 1752 2010-04-22 09:11:05Z dmessina $

use Modern::Perl;
use Bio::DB::EUtilities;
use Getopt::Long;

my ($besthit);
GetOptions( 'besthit' => \$besthit );

my $usage = "acc2date.pl - get the last-modified date for NCBI IDs

Usage: acc2date.pl <file_of_IDs>

where file_of_IDs has one ID per line.

Options:
--besthit    input file is output from besthit, which looks like:

cluster0	all_c751	gi|136358448|gb|EBN62812.1|	86.64	247	26	2	1	720	45	291	2e-88 326

With --besthit, output will be the besthit line with the hit's date appended.
";
die $usage unless scalar @ARGV == 1;

# slurp lines from the file
open( my $fh, $ARGV[0] ) or die "couldn't open $ARGV[0]: $!\n";
my @lines = <$fh>;
chomp @lines;
my @ids;

my %id_map;
my %besthits;
if ($besthit) {
    foreach my $line (@lines) {
        my @fields = split "\t", $line;

        if ($fields[1] eq '- none -') {
            $id_map{ $fields[0] } = { 'date' => '',
                                      'line' => $line,};
            next;
        }
        
        die ("This input line doesn't look right:\n", $line, "\n")
          if scalar @fields != 13;

        my $raw_id_field = $fields[2];
        my $id;
        if ($raw_id_field =~ /gi\|(\w+?)\|/) {
            $id = $1;
            push @ids, $id;
        }
        else {
            die "couldn't get gi from $raw_id_field\n";
        }

        # use ID as key in besthit lookup
        push @{$besthits{$id}}, $line;
    }
}
else {
    @ids = @lines;
}


my $eutil = Bio::DB::EUtilities->new(
    -eutil  => 'esummary',
    -db     => 'protein',
    -email  => 'david.messina@sbc.su.se',
    -retmax => 2500,
    -id     => \@ids,
);


while ( my $ds = $eutil->next_DocSum ) {
    my $id = $ds->get_id;
    my ($cdate) = $ds->get_contents_by_name('UpdateDate');
    die "no date for $id!" if !defined($cdate);

    if ($besthit) {
        # lookup the line that goes with the ID
        # my @keys = grep /$id/, keys %besthits;
        # die join("\t", "multiple matches for $id", @keys), "\n" if scalar @keys > 1;
        # my $line = $besthits{$keys[0]};

        my $linesref = $besthits{$id};
        my @lines = @$linesref;
        die "no lines for $id!" if !defined($lines[0]);

        foreach my $line (@lines) {
            my ($cluster) = split("\t", $line, 2);
            
            # save the stuff to the ID map
            if (exists $id_map{ $cluster }) {
                warn "$cluster already has an entry!\n";
            }
            else {
                $id_map{ $cluster } = { 'date' => $cdate,
                                        'line' => $line,};
            }
        }
    }
    else {
        $id_map{ $id } = $cdate;
    }
}



if ($besthit) {
    say join( "\t", $id_map{$_}->{'line'}, $id_map{$_}->{'date'} ) for sort by_num keys %id_map;
}
else {
    say join( "\t", $_, $id_map{$_} ) for sort by_num keys %id_map;    
}


sub by_num {
    ( $a =~ /(\d+)/ )[0] <=> ( $b =~ /(\d+)/ )[0]
}
