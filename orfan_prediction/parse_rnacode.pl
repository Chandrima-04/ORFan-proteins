#!/usr/bin/perl
# $Id: parse_rnacode.pl 3117 2011-09-15 09:56:51Z dmessina $
# created by Dave Messina on 2010-05-24

use Modern::Perl;
use Getopt::Long;
use File::Basename;
use Cluster;
use ClusterGroup;

my $pval = 0.05;
my ($cl, $noheaders);
my $len = 0;
GetOptions('pval:f'    => \$pval,
           'cl:s'      => \$cl,
           'len:i'     => \$len,
           'noheaders' => \$noheaders,
           );

my @outlines;

my $usage = "parse_rnacode.pl - parses output from RNAcode, grabs best hit
                   from each file, and prints all results with
                   P-value <= the cutoff
                   
Usage: parse_rnacode.pl --cl <.cl file> --pval 0.15 cluster1.rnacode [cluster2.rnacode] ...

where .cl file is a ClusterGroup file

Options:
--pval n     set the P-value cutoff (default 0.05)
--len  n     set the HSSP length cutoff (default 0)
--noheaders   suppress column headers in output
";
die $usage unless @ARGV;

if (! -e $cl) {
    say STDERR "\nERROR! could not find $cl: $!", "\n";
    die $usage;
}

my $cluster_group = ClusterGroup->new('-load' => $cl);

foreach my $file (@ARGV) {
    open (my $fh, '<', $file) or die "couldn't open $file: $!\n";
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );

    LINE:while (my $line = <$fh>) {
        chomp $line;
        
        $line =~ s/^\s+//; #strip leading whitespace
        my @fields = split /\s+/, $line;
        if ( scalar @fields == 11 && $fields[0] == '0' ) {

            # does it pass our pval cutoff?
            next LINE if ($fields[10] > $pval);
            
            # does it pass out HSSP length cutoff?
            next LINE if ($fields[3] < $len);
            
            # lookup number of seqs
            my $cluster_size;
            if ($base =~ /cluster(\d+)/) {
                my $cluster_id = $1;
                my $cluster = $cluster_group->get_cluster('-id' => $cluster_id);
                $cluster_size = $cluster->size;
            }
            else {
                die "couldn't get cluster id for $base\n";
            }

            my $outline = join("\t", "$base   ", $cluster_size, $line);
            push @outlines, $outline;                
        }
    }

}

unless ($noheaders) {
    say '                                #seqs   #       Strand  Frame  Length   From    To      Name                    Start   End       Score      P';
    say '                                =================================================================================================================';
}

foreach my $result (sort by_pval @outlines) {
    say STDOUT $result;
}

sub by_pval {
    my @a_fields = split(/\s+/, $a);
    my @b_fields = split(/\s+/, $b);

    $a_fields[12] <=> $b_fields[12]
        ||
    $b_fields[11] <=> $a_fields[11]
}
