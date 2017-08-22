#!/usr/bin/perl
# $Id: add_seg_stats.pl 3123 2011-09-16 09:35:29Z dmessina $
# created by Dave Messina on 2011-07-06

use Modern::Perl;
use Getopt::Long;
use File::Basename;
use File::Spec;
use Bio::SeqIO;
use Bio::AlignIO;
use List::Util qw(min max sum);
use POSIX;

my ($segdir, $aa_alndir, $nt_alndir, $rnacode_dir, $debug);
GetOptions('segdir:s'      => \$segdir,
           'aa_alndir:s'   => \$aa_alndir,
           'nt_alndir:s'   => \$nt_alndir,
           'rnacode_dir:s' => \$rnacode_dir,
           'debug'         => \$debug,);

my $usage = "add_seg_stats.pl - add seg stats to rnacode summary file

Usage: add_seg_stats.pl rnacode_summary.txt --segdir path/to/segfiles --aa_alndir path/to/aln
        --nt_alndir path/to/aln
where:
aa_alndir contains clustal protein alignments, one for each cluster
nt_alndir contains clustal nucleotide alignments, one for each cluster
rnacode_dir contains eps directories for each cluster named like cluster1.fa_eps

NOTE: aligments must be in clustal (.aln) format, and the filenames must end
      in '.aln.txt' (nt) and '_aa.aln.txt'.
";
die $usage unless @ARGV && -d $segdir && -d $aa_alndir && -d $nt_alndir
    && -d $rnacode_dir;


foreach my $file (@ARGV) {
    open (my $fh, '<', $file) or die "couldn't open $file: $!\n";
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my @unsorted_output;

    while (my $line = <$fh>) {
        chomp $line;

        # skip header lines
        next unless ($line =~ /^\w+/);
        
        my @fields = split /\s+/, $line;
        
        # validate line format
        unless ( scalar @fields == 13 ) {
            die "wrong number of fields on line:\n$line\n";
        }
        
        # validate clusterID in first field
        my ($fileroot, $cluster_id, $block);
        if ($fields[0] =~ /(cluster(\d+).+?)\_/) {
            ($fileroot, $cluster_id) = ($1, $2);
        }
        else { die "couldn't get cluster id from $fields[0]\n"; }

        my $seq_count = $fields[1];
        
        # fetch and parse seg files
        my @seghi_stats = parse_seg($cluster_id, '.fa_aa_seghi.fa', $seq_count);
        my @seglo_stats = parse_seg($cluster_id, '.fa_aa_seglo.fa', $seq_count);
        my @seg_stats;
        
        # take whichever's longer
        if (!exists $seglo_stats[0] || ($seghi_stats[5] >= $seglo_stats[5])) {
            @seg_stats = ('hi', @seghi_stats);
        }
        else { @seg_stats = ('lo', @seglo_stats); }

        # get alignment (ORF) length
        my $aln_length = aln_length($aa_alndir, $cluster_id);
        my $log_aln_length = sprintf("%.2f", log10($aln_length));

        # log10 of RNAcode score
        my $log_rnacode = sprintf("%.2f", log10($fields[11]));
            
        # get rid of extraneous columns
        my @new_fields = ($fields[0], $fields[1], $fields[11], $log_rnacode, 
            $fields[12]);

        # ad hoc combined score - $seg_stats[2] is avg complexity
        my $combo_score = $log_aln_length + $log_rnacode + $seg_stats[2];

        # add links
        my $nt_aln_path  = File::Spec->catfile($nt_alndir, "${fileroot}.aln.txt");
        my $aa_aln_path  = File::Spec->catfile($aa_alndir, "${fileroot}_aa.aln.txt");
        my $eps_path     = File::Spec->catfile($rnacode_dir, "${fileroot}_eps", 'hss-0.eps');
        my $nt_aln_link  = '<a href="'. $nt_aln_path . '">DNA</a>';
        my $aa_aln_link  = '<a href="'. $aa_aln_path . '">AA</a>';
        my $eps_link     = '<a href="'. $eps_path . '">EPS</a>';

        # collect all the fields in a new array and store it
        my @output = (@new_fields, $aln_length, $log_aln_length, @seg_stats,
                      $combo_score, $nt_aln_link, $aa_aln_link, $eps_link);                      
        push @unsorted_output, \@output;
    }
    
    # write markdown metadata and table header
    md_metadata($file);
    table_header();
    
    # sort by combo score and write out
    my @sorted_output = sort { $b->[16] <=> $a->[16] } @unsorted_output;
    foreach my $row (@sorted_output) {
        say join "\t|", @$row;
    }
}

sub parse_seg {
    my ($cluster_id, $suffix, $seq_count) = @_;
    my (@ids, @complexity, @length);
    my %seg_data;

    # open seg file
    my $segfile = 'cluster' . $cluster_id . $suffix;
    my $segpath = File::Spec->catfile($segdir, $segfile);
    my $seg_io  = Bio::SeqIO->new(-format => 'fasta',
                                  -file   => $segpath);

    # grab seg and length data for each seq in the cluster
    while(my $seq_obj = $seg_io->next_seq) {
        my $raw_id   = $seq_obj->primary_id;
        my ($id, $start, $end);
        if ($raw_id =~ /(\w+)\((\d+)\-(\d+)\)/) {
            ($id, $start, $end) = ($1, $2, $3);
        }
        my $length = $end - $start + 1;

        my $raw_desc = $seq_obj->desc;
        my $complexity;
        if ($raw_desc =~ /complexity=(.+)\s\(/) {
            $complexity = $1;
        }
        # say join "\t", $cluster_id, $id, $start, $end, $length, $complexity;

        if (exists $seg_data{$id} && $seg_data{$id}->{'length'} > $length) {
            my $lendiff  = $seg_data{$id}->{'length'} - $length;
            my $cmplxdiff = $seg_data{$id}->{'complexity'} - $complexity;
            if ($length > 1) { say STDERR "len > by $lendiff, cmplx diff by $cmplxdiff"; }
            next;
        }
        else {
            $seg_data{$id}->{'complexity'} = $complexity;
            $seg_data{$id}->{'length'}     = $length;
        }
    }

    # confirm I got data on all the sequences
    unless (scalar keys %seg_data == $seq_count) {
        warn "didn't get data on all the seqs in cluster $cluster_id\n";
    }

    foreach my $id (keys %seg_data) {
        push @complexity, $seg_data{$id}->{'complexity'};
        push @length, $seg_data{$id}->{'length'};
    }

    # check that there's actually data to summarize
    if (scalar @complexity == 0) {
        return;
    }

    my $min_complexity = min @complexity;
    my $max_complexity = max @complexity;
    my $avg_complexity = sprintf("%.2f", (sum @complexity) / scalar @complexity);
    my $maxcmplx_count = grep($_ == $max_complexity, @complexity);
    my $min_length = min @length;
    my $max_length = max @length;
    my $avg_length = sprintf("%.2f", (sum @length) / scalar @length);
    my $maxlen_count = grep($_ == $max_length, @length);
    my @seg_stats = ($min_complexity, $avg_complexity, $max_complexity,
                    $maxcmplx_count, $min_length, $avg_length,
                    $max_length, $maxlen_count);
    return @seg_stats;
}

sub aln_length {
    my ($aa_alndir, $cluster_id) = @_;
    my $file = 'cluster' . $cluster_id . '.fa_aa.aln.txt';
    my $path = File::Spec->catfile($aa_alndir, $file);
    my $aln_io = Bio::AlignIO->new(-file => $path,
                                 -format => 'clustalw');
    my $aln = $aln_io->next_aln;
    my $length = $aln->length;
    return $length;
}

# spit out metadata for Markdown to make HTML from
sub md_metadata {
    my ($file) = @_;
    say 'Title: ', "Summary report for $file";
    say 'CSS: ', 'summaryreport.css';
    say 'Author: ', 'Dave Messina', "\n\n";
}

sub table_header {
    # overcolumns group columns by type
    my $empty = ' ';
    my @overcolumns = ($empty, $empty, 'RNAcode ||',
                       $empty, $empty,
                       'Seg complexity ||||',
                       'Seg length |||',
                       $empty,
                       $empty, $empty, $empty,);
    my @columns = ('Cluster', '# seqs', 'Score', 'log RNAcode score', 'P',
                  'ORF length (aa)', 'log ORF length',
                  'complexity', 'min', 'avg', 'max', '# seqs with max complexity',
                  'min', 'avg', 'max', '# seqs with max length',
                  'composite score',
                  'DNA alignment', 'protein alignment', 'RNAcode figure',);
    my $centered = ':-:';
    my $right    = '-';
    my @dividers = ($centered, $centered, $centered, $centered, $centered,
                    $centered, $centered,
                    $centered, $centered, $centered, $centered, $centered,
                    $centered, $centered, $centered, $centered,
                    $centered,
                    $centered, $centered, $centered,);
    # my @dividers = ($centered, $right, $right, $right, $right,
    #                 $right, $right,
    #                 $centered, $right, $right, $right, $right,
    #                 $right, $right, $right, $right,
    #                 $right,
    #                 $centered, $centered, $centered,);
    say '|', join('|', @overcolumns), '|';
    say '|', join('|', @columns    ), '|';
    say '|', join('|', @dividers   ), '|';
}