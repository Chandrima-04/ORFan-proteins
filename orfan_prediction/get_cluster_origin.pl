#!/usr/bin/perl
# $Id$
# created by Dave Messina on 2011-08-22

use Modern::Perl;
use Bio::SeqIO;
use Statistics::Descriptive;

my $usage = "get_cluster_origin.pl - determine which samples the cluster's seqs are composed of
Usage: get_cluster_origin.pl ~/Projects/virus/results/20110706/ten_summary_sortedbylen+cmplx_cluster_names.txt
";
die $usage unless @ARGV;

my $fasta_dir = '/Users/dave/Projects/virus/results/20110706/ten_prot_fastas';

my %seen_clusters;

# print header
print join("\t", 'cluster',
    'CSF count', 'CSF %',
    'contig count', 'contig %',
    'feces count', 'feces %',
    'mucus count', 'mucus %',
    'serum count', 'serum %',
    'tonsil count', 'tonsil %',
    'DNA %', 'RNA %','unk %'), "\n";

while (my $line = <>) {

    # grab cluster name and convert to a fasta filename
    chomp $line;
    my $raw_cluster = $line;
    my ($cluster_id, $rest) = split /\./, $raw_cluster;
    my $cluster_filebase = $cluster_id . '_aa.fa';
    my $fasta_file = $fasta_dir . '/' . $cluster_filebase;

    # these lines are a list of blocks, so the same cluster can appear
    # more than once. check for that.
    if ($seen_clusters{$cluster_id}) { next; }
    else { $seen_clusters{$cluster_id}++; }
    
    # open fasta file and grab seq IDs
    my %sources = ('mucus' => 0,
                   'serum' => 0,
                   'cerebrospinal fluid' => 0,
                   'feces' => 0,
                   'tonsil' => 0,
                   'contig' => 0);
    my %types = ('DNA' => 0,
                'RNA' => 0,
                'unk' => 0,);
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
                                -file   => $fasta_file);
    my $seq_count;
    while (my $entry = $seqio->next_seq) {
        $seq_count++;
        my ($source, $type) = id2source( $entry->display_id() );
        $sources{$source}++;
        $types{$type}++;
    }
    
    # print summary
    print $cluster_id, "\t";
    foreach my $key (sort keys %sources) {
        my $source_pc = sprintf("%.2f", ($sources{$key}/$seq_count)*100);
        print join("\t", $sources{$key}, $source_pc), "\t";
    }
    my $dna_pc = sprintf("%.2f", ($types{'DNA'}/$seq_count)*100);
    my $rna_pc = sprintf("%.2f", ($types{'RNA'}/$seq_count)*100);
    my $unk_pc = sprintf("%.2f", ($types{'unk'}/$seq_count)*100);
    print join("\t", $dna_pc, $rna_pc, $unk_pc), "\n";
}

sub id2source {
    my ($id) = @_;
    
    my $id_root = substr($id, 0, 9);
    
    my %sources = (
        'ER8QEOW01'              =>   'mucus',
        'ER8QEOW02'              =>   'mucus',
        'FC8LRL301'              =>   'mucus',
        'FC8LRL302'              =>   'mucus',
        'FCPU0RF01'              =>   'serum',
        'FCPU0RF02'              =>   'serum',
        'FPFNBIF01'              =>   'serum',
        'FPFNBIF02'              =>   'serum',
        'FS22EC101'              =>   'serum',
        'FS22EC102'              =>   'serum',
        'FSTRRC101'              =>   'serum',
        'FTSPZO101'              =>   'cerebrospinal fluid',
        'FTSPZO102'              =>   'cerebrospinal fluid',
        'GB3LKKR01'              =>   'feces',
        'GB3LKKR02'              =>   'feces',
        'GB7HT0Z01'              =>   'tonsil',
        'GB7HT0Z02'              =>   'tonsil',
    );
    
    # determine molecule type from last digit of source ID
    my ($source, $type);
    if ($sources{$id_root}) {
        $source = $sources{$id_root};
        my $type_digit = substr($id_root, 8, 1);
        $type = ($type_digit == 1) ? 'DNA' : 'RNA';

        # NOTE! does not apply to tonsil -- both libs are RNA!
        if ($source eq 'tonsil') {
            $type = 'RNA';
        }
    }
    else {
        $source = 'contig';
        $type   = 'unk';
    }

    return ($source, $type); 
}