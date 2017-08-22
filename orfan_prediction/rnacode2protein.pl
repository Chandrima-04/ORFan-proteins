#!/usr/bin/perl
# $Id: rnacode2protein.pl 3115 2011-09-15 09:28:15Z dmessina $
# created by Dave Messina on 2011-06-23

use Modern::Perl;
use Getopt::Long;
use File::Basename;
use File::Spec;
use Bio::AlignIO;

my ($cl, $alndir, $outdir, $debug, $refonly, $refgaps);
GetOptions('cl:s'     => \$cl,
           'alndir:s' => \$alndir,
           'outdir:s' => \$outdir,
           'debug'    => \$debug,
           'refonly'  => \$refonly,);

my $usage = "rnacode2protein.pl - extract predictions protein sequences from RNAcode output

Usage: rnacode2protein.pl rnacode_summary.txt --cl my.cl --alndir path/to/aln --outdir path/to/outdir

where .cl file is a ClusterGroup file
and alndir points to a directory containing the ClustalW format multiple alignments for each cluster

--refonly only write out the translation of the reference sequence (first seq in the alignment)
--debug turn on some debugging output

";
die $usage unless @ARGV && -e $cl && -d $alndir && -d $outdir;


foreach my $file (@ARGV) {
    open (my $fh, '<', $file) or die "couldn't open $file: $!\n";
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );

    while (my $line = <$fh>) {
        chomp $line;

        # skip header lines
        next unless $line =~ /^\w+/;
        
        my @fields = split /\s+/, $line;
        
        # validate line format
        unless ( scalar @fields == 13 ) {
            die "wrong number of fields on line:\n$line\n";
        }
        
        # validate clusterID in first field
        my ($fileroot, $cluster_id);
        if ($fields[0] =~ /(cluster(\d+).+?)\_/) {
            ($fileroot, $cluster_id) = ($1, $2);
        }
        else { die "couldn't get cluster id from $fields[0]\n"; }

        # open alignment file
        my $alnfile = $fileroot . '.aln';
        my $alnpath = File::Spec->catfile($alndir, $alnfile);
        my $aln_io = Bio::AlignIO->new(
            '-file'   => $alnpath,
            '-format' => 'clustalw',
            '-displayname_flat' => 1);
        my $aln = $aln_io->next_aln;
        unless (defined $aln) {
            die "couldn't get alignment from $alnpath\n";
        }

        # info to the user
        say STDERR "cluster $cluster_id";

        # extract refseq name
        my $refseq_name = $fields[8];
        $refseq_name =~ s/\/\d+\-\d+$//;

################ REFSEQ ONLY ###################
        # extract refseq's DNA
        my $refseq_dna = $aln->get_seq_by_id($refseq_name);

        # pull ref seq and adjust for frame
        my $dna_string = $refseq_dna->seq;
        $dna_string =~ s/-//g; # remove gaps
        my $frame = $fields[4] - 1;
        my $is_rev = $fields[3] eq '-' ? 1 : 0;

        # make new object, revcomp if - strand
        my $gapless_refseq_dna_obj = $refseq_dna->clone(-seq => $dna_string);
        if ($is_rev) {
            $gapless_refseq_dna_obj = $gapless_refseq_dna_obj->revcom();
        }

        # reset coords and translate in correct frame
        $gapless_refseq_dna_obj->start(1);
        $gapless_refseq_dna_obj->end($gapless_refseq_dna_obj->length);
        my $refseq_aa_obj = $gapless_refseq_dna_obj->translate(-frame => $frame);

        # extract just the HSS part of the sequence
        my $hss_protein_string = $refseq_aa_obj->subseq($fields[6], $fields[7]);

        # extract the whole ORF containing the HSS
        my ($protein_orf_obj, $prot_orf_start, $prot_orf_end) =
            whole_orf_from_coords($refseq_aa_obj, $fields[6], $fields[7]);
        my $orf_protein_string = $protein_orf_obj->seq;

        # make prot aln from orf to map back to DNA coords of refseq
        my ($dna_start, $dna_end) = map_prot_to_dna($refseq_aa_obj,
             $refseq_dna, $prot_orf_start,
             $prot_orf_end, $frame, $is_rev);
        if ($debug) {
            say join "\t", "aa_start:$prot_orf_start",
                           "dna_start:$dna_start",
                           "aa_end:$prot_orf_end",
                           "dna_end:$dna_end";
        }
        
        # write out protein seq
        if ($debug) {
            say STDOUT '>', 'cl', $cluster_id, '_', $refseq_aa_obj->display_id;
            say $orf_protein_string;
        }

        # make orf-only dna object
        my $ref_orf_dna_obj = $refseq_dna->trunc($dna_start, $dna_end);
        if ($is_rev) {
            $ref_orf_dna_obj = $ref_orf_dna_obj->revcom();
        }
        my $ref_orf_dna_string = $ref_orf_dna_obj->seq;

        # write out DNA string (debug only)
        if ($debug) {
            my $orf_dna_string = $refseq_dna->subseq($dna_start, $dna_end);
            say STDERR '>', 'cl', $cluster_id, '_', $refseq_name,
                ' with gaps', "\n", $orf_dna_string;
        }

        # stop here if only doing refseq
        next if $refonly;
        
################ REFSEQ ONLY ###################

        # slice out the orf part of the alignment
        # 3rd arg '1' is boolean to preserve gap-only columns
        my $whole_orf_aln = $aln->slice($dna_start, $dna_end, 1);

        # reverse comp whole alignment if necessary
        if ($is_rev) {
            $whole_orf_aln = revcom_aln($whole_orf_aln);
        }

        # remove gap columns and translate one codon at a time
        my @ref_orf = split '', $ref_orf_dna_string;
        my (@gap_cols, @coords, @codon, %protseqs);
        my $pos_in_this_codon = 0;
        for (my $i = 0; $i < scalar @ref_orf; $i++) {
            my $nt = $ref_orf[$i];

            # mark gaps
            if ($nt eq '-') {
                push @gap_cols, $pos_in_this_codon;
            }
            else {
                push @codon, $nt;
                push @coords, $i;
            }
            $pos_in_this_codon++;

            # 3 non-gap cols is a codon, so slice it out and translate
            if (scalar @coords == 3) {
                
                # slice out that codon
                my $codon_start = $coords[0]+1; # slice is 1-based
                my $codon_end   = $coords[2]+1;
                my $codon_aln = $whole_orf_aln->slice($codon_start, $codon_end, 1);
                
                # remove the gap columns ( while shift gap_cols)
                while (my $pos = shift @gap_cols) {
                    my $new_codon_aln = $codon_aln->remove_columns([$pos, $pos]);
                    $codon_aln = $new_codon_aln;
                }
                
                # translate codon in every seq
                for (my $i = 1; $i <= $codon_aln->num_sequences; $i++) {
                    my $codon = $codon_aln->get_seq_by_pos($i);
                    my $aa;
                    # all gaps becomes a gap
                    if ($codon =~ /\-{3}/)      { $aa = '-'; }
                    # 1 or 2 gaps becomes an unknown aa
                    elsif ($codon =~ /\-{1,2}/) { $aa = 'X'; }
                    # 0 gaps becomes a known aa
                    else { $aa = $codon->translate(-frame => 0); }
                    # add to the growing protein sequence
                    
                    $protseqs{$aa->display_id}->{'seq'} .= $aa->seq;
                    $protseqs{$aa->display_id}->{'pos'}  = $i;
                }

                # cleanup
                undef @gap_cols;
                undef @codon;
                undef @coords;
                $pos_in_this_codon = 0;
            }
        }
        
        # new prot alignment, same seq order as dna alignment
        my $prot_align = Bio::SimpleAlign->new();
        foreach my $seq (sort { $protseqs{$a}{'pos'} <=> $protseqs{$b}{'pos'} } keys %protseqs) {
            my $new_aa_obj = Bio::LocatableSeq->new(-display_id => $seq,
                -alphabet => 'protein',
                -strand   => 1,
                -seq      => $protseqs{$seq}{'seq'}, );
            $prot_align->add_seq($new_aa_obj);
        }
        
        # write out
        my $outfile = $fileroot . '_aa.aln';
        my $outpath = File::Spec->catfile($outdir, $outfile);
        my $out = Bio::AlignIO->new(-file             => ">$outpath",
                                    -format           => 'clustalw',
                                    -displayname_flat => 1);
        $out->write_aln($prot_align);
    }    
}

# note: this strips off all the '*' stop codons, from the beginning and
# the end of the ORF
sub whole_orf_from_coords {
    my ($seq_obj, $start, $end) = @_;
    my $raw_seq = $seq_obj->seq;
    my @raw_chars = split '', $raw_seq;
    my ($orf_start, $orf_end);
    
    # extend to end
    my $end_char = $seq_obj->subseq($end, $end);
    if ($end_char eq '*') {
        $orf_end = $end - 1;
    }
    else {
        for (my $i = $end; $i < @raw_chars; $i++) {
            if ($raw_chars[$i] eq '*') {
                $orf_end = $i;
                last;
            }
        }
        if (!defined $orf_end) {
            # if no stop char, end of orf is the end of the sequence
            $orf_end = scalar @raw_chars;
        }
    }
    
    # extend to start
    my $start_char = $seq_obj->subseq($start, $start);
    if ($start_char eq '*') {
        $orf_start = $start + 1;
    }
    else {
        for (my $i = $start-2; $i >= 0; $i--) {
            if ($raw_chars[$i] eq '*') {
                $orf_start = $i + 2;
                last;
            }
        }
        if (!defined $orf_start) {
            # if no stop char, start of orf is the start of the sequence
            $orf_start = 1;
        }
    }

    # use coordinates to extract the whole orf seq
    my $protein_orf_obj = $seq_obj->trunc($orf_start, $orf_end);
    return ($protein_orf_obj, $orf_start, $orf_end);
}

sub map_prot_to_dna {
    my ($aa_nogap_obj, $dna_obj, $aa_start, $aa_end, $frame, $is_rev) = @_;
    my ($dna_start, $dna_end);

    my $aa_nogap_len = $aa_nogap_obj->length;
    my $aa_nogap_string = $aa_nogap_obj->seq;

    my $dna_string;
    if ($is_rev) { $dna_string = $dna_obj->revcom->seq; }
    else         { $dna_string = $dna_obj->seq; }

    my $dna_len = length($dna_string);

    my $start_offset = 0;

    my $dna_pos = 1; # 1-based, DNA coordinates
    
    # walk along the aa seq one aa at a time and follow on the dna seq
    for (my $i = 0; $i < $aa_nogap_len; $i++) {
        my $char = substr($aa_nogap_string, $i, 1);
        
        # 3 times, once for each position in the codon
        my (@codon, @dna_coords);
        for (1..3) {
            my ($nt, $pos) = get_next_char($dna_string, $dna_pos);
            push @codon, $nt;
            push @dna_coords, $pos;
            $dna_pos = $pos + 1; # move forward on position to get next char
        }
        
        my $codon_string  = join '', @codon;
        my $coords_string = join ',', @dna_coords;
        my $aa_pos = $i+1;

        say join "\t", $aa_pos, $char, $codon_string, $coords_string if $debug;

        if ($is_rev) {
            $codon_string = join '', reverse(@codon);
            $codon_string =~ tr/ATGCNatgcn/TACGNtacgn/;
        }
        
        if ($aa_pos == $aa_start) {
            $dna_start = $dna_coords[0];
        }
        if ($aa_pos == $aa_end) {
            $dna_end = $dna_coords[2];
        }
        if (defined $dna_start && defined $dna_end) {
            
            if ($is_rev) {
                $dna_start = $dna_len - $dna_start + 1;
                $dna_end   = $dna_len - $dna_end   + 1;
                
                # exchange start and end since for slice start must < end
                if ($dna_start > $dna_end) {
                    ($dna_start, $dna_end) = ($dna_end, $dna_start);
                }
            }

            # adjust DNA pos based on frame and strand
            my $adjustment = $frame;
            if ($is_rev) { $adjustment = '-' . $adjustment; }
            $dna_start += $adjustment;
            $dna_end   += $adjustment;
            
            # if start is on a gap, add one
            # -1 because substr takes offset not abs pos
            # must reverse dna string before this part because dna_start
            # and dna_end positions are relative to the strand
            if ($is_rev) { $dna_string = reverse($dna_string); }
            my $first_char = substr( $dna_string, $dna_start-1, 1);
            if ($first_char eq '-') { $dna_start++; }

            return ($dna_start, $dna_end);
        }
    }


}

sub get_next_char {
    my ($seq, $cur_pos) = @_;
    
    # cur_pos is in DNA, 1-based coordinates
    # $i is in array, 0-based coordinates
    
    my @seq = split '', $seq;
    
    for (my $i = $cur_pos - 1; $i < scalar @seq; $i++) {
        my $nt = $seq[$i];
        if ($nt eq '-') {
            next;
        }
        else {
            my $new_pos = $i + 1; # convert back to DNA-based coordinates
            return ($nt, $new_pos);
        }
    }
}

sub revcom_aln {
    my ($whole_orf_aln) = @_;
    
    for (my $i = 1; $i <= $whole_orf_aln->num_sequences; $i++) {
        my $seq_obj = $whole_orf_aln->get_seq_by_pos($i);
        my $rev_seq_obj = $seq_obj->revcom();
        $whole_orf_aln->remove_seq($seq_obj);
        $whole_orf_aln->add_seq(-SEQ => $rev_seq_obj, -ORDER => $i);
    }

    return $whole_orf_aln;
}