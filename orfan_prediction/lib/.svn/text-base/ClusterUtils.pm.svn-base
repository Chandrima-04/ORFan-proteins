package ClusterUtils;

use strict;
use warnings;
use File::Temp ();

use base 'Exporter';
our @EXPORT = qw(fastacmd blastdbcmd sfetch normalize_seqid);

# ClusterUtils - convenience functions for creating Cluster objects

# fastacmd - fetches seqs from a formatdb database and returns a
# filehandle to those seqs in FASTA format
#
# Arguments:
# - the path to a formatdb database
# - a reference to an array of sequence IDs
#
# Returns:
# a filehandle which can be passed to Bio::SeqIO
sub fastacmd {
    my ($seqdb, $seqarrayref) = @_;
    
    # get seqs using fastacmd
    my $cluster_size = scalar @{ $seqarrayref };
    my @args;
    if ($cluster_size < 1000) {
        my $id_string = join(',', @{ $seqarrayref });
        @args = ('fastacmd', '-d', $seqdb, '-s', $id_string);
    }
    # put seqids in tempfile for big clusters
    else {
        # tempfile will not be auto-deleted, but it's usually written to
        # /tmp or similar, so the OS will get rid of it eventually.
        my $tempfh = File::Temp->new('UNLINK' => 0);
        print $tempfh join("\n", @{ $seqarrayref }), "\n";
        @args = ('fastacmd', '-d', $seqdb, '-i', $tempfh->filename);
    }
    
    my $pid = open(my $seqfh, '-|', @args);
    if (!defined($pid)) { die "can't fastacmd: $!"; }

    
    return $seqfh;
}


# blastdbcmd - fetches seqs from a makeblastdb (BLAST+) database and returns a
# filehandle to those seqs in FASTA format
#
# Arguments:
# - the path to a makeformatdb database
# - a reference to an array of sequence IDs
#
# Returns:
# a filehandle which can be passed to Bio::SeqIO
sub blastdbcmd {
    my ($seqdb, $seqarrayref) = @_;
    
    # get seqs using fastacmd
    my $cluster_size = scalar @{ $seqarrayref };
    my @args;
    if ($cluster_size < 1000) {
        my $id_string = join(',', @{ $seqarrayref });
        @args = ('blastdbcmd', '-db', $seqdb, '-entry', $id_string);
    }
    # put seqids in tempfile for big clusters
    else {
        # tempfile will not be auto-deleted, but it's usually written to
        # /tmp or similar, so the OS will get rid of it eventually.
        my $tempfh = File::Temp->new('UNLINK' => 0);
        print $tempfh join("\n", @{ $seqarrayref }), "\n";
        @args = ('blastdbcmd', '-db', $seqdb, '-entry_batch', $tempfh->filename);
    }
    
    my $pid = open(my $seqfh, '-|', @args);
    if (!defined($pid)) { die "can't blastdbcmd: $!"; }

    
    return $seqfh;
}

# sfetch - fetches seqs from a esl-sfetch database and returns a
# filehandle to those seqs in FASTA format
#
# Arguments:
# - the path to an sfetch database
# - a reference to an array of sequence IDs
#
# Returns:
# a filehandle which can be passed to Bio::SeqIO
sub sfetch {
    my ($seqdb, $seqarrayref) = @_;
    
    # get seqs using fastacmd
    my $cluster_size = scalar @{ $seqarrayref };
    my @args;

    # tempfile will not be auto-deleted, but it's usually written to
    # /tmp or similar, so the OS will get rid of it eventually.
    my $tempfh = File::Temp->new('UNLINK' => 0);
    print $tempfh join("\n", @{ $seqarrayref }), "\n";
    @args = ('esl-sfetch', '-f', $seqdb, $tempfh->filename);
    
    my $pid = open(my $seqfh, '-|', @args);
    if (!defined($pid)) { die "can't sfetch: $!"; }

    
    return $seqfh;
}


# convenience method to fix seqIDs
# - really this should be in Bio::SeqUtils
# if seqID is lcl|A12345 convert it to simply A12345
# might want in the future to detect and convert other common
# seqID prefixes like gb| 
sub normalize_seqid {
    my ($seq) = @_;
    
    my $current_id = $seq->display_id();
    if ( $current_id =~ /^lcl\|(.+)/ ) {
        $seq->display_id($1);
    }
    
    return $seq;
}

1;