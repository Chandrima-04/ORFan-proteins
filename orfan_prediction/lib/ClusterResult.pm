package ClusterResult;

# $Id: ClusterResult.pm 2940 2011-04-19 14:16:48Z dmessina $
# created by Dave Messina on 2011-04-18

use strict;
use warnings;

sub new {
    my ( $class, %args ) = @_;

    my $object = {
        'SS'   => 0,
        'SD'   => 0,
        'DS'   => 0,
        'DD'   => 0,
        'corrected' => 0,
        'DS_seqs' => {},
        'SD_seqs' => {},
    };
    bless( $object, $class );

    return $object;
}

sub SS {
    my ( $self, $value ) = @_;
    if ( defined $value ) { $self->{'SS'} += $value; }
    else                  { return $self->{'SS'}; }
}

sub SD {
    my ( $self, $value ) = @_;
    if ( defined $value ) { $self->{'SD'} += $value; }
    else                  { return $self->{'SD'}; }
}

sub DS {
    my ( $self, $value ) = @_;
    if ( defined $value ) { $self->{'DS'} += $value; }
    else                  { return $self->{'DS'}; }
}

sub DD {
    my ( $self, $value ) = @_;
    if ( defined $value ) { $self->{'DD'} += $value; }
    else                  { return $self->{'DD'}; }
}

sub jaccard {
    my ($self) = @_;

    my $jaccard = $self->SS / ( $self->SS + $self->DS + $self->SD );
    $jaccard = sprintf( "%.2f", $jaccard );

    return $jaccard;
}

sub rand {
    my ($self) = @_;

    my $total = $self->SS + $self->SD + $self->DS + $self->DD;
    my $rand  = ( $self->SS + $self->DD ) / $total;
    $rand = sprintf( "%.2f", $rand );

    return $rand;
}

sub folkes_mallows {
    my ($self) = @_;

    my $fmi =
      $self->SS / sqrt( ( $self->SS + $self->DS ) * ( $self->SS + $self->SD ) );
    $fmi = sprintf( "%.2f", $fmi );

    return $fmi;
}

sub in_cluster {
    my ($self) = @_;

    my $in_cluster = $self->SS / ( $self->SS + $self->SD );
    $in_cluster = sprintf( "%.2f", $in_cluster );

    return $in_cluster;
}

sub seq_diff {
    my ($self) = @_;

    my $seq_diff = $self->DD / ( $self->DS + $self->DD );
    $seq_diff = sprintf( "%.2f", $seq_diff );

    return $seq_diff;
}

sub cluster_start {
    my ( $self, $value ) = @_;
    if ( defined $value ) { $self->{'start'} = $value; }
    else                  { return $self->{'start'}; }
}

sub cluster_end {
    my ( $self, $value ) = @_;
    if ( defined $value ) { $self->{'end'} = $value; }
    else                  { return $self->{'end'}; }
}

sub to_string {
    my ($self) = @_;
    
    my $string = join("\t", $self->cluster_start, $self->cluster_end,
    $self->SS, $self->SD, $self->DD, $self->DS);
    
    return $string;
}

# halve all the raw stats since I use a full matrix
sub correct {
    my ($self) = @_;

    foreach my $method ('SS', 'SD', 'DS', 'DD') {
        my $value = $self->$method;
        $self->{$method} = $value/2;
    }
    $self->{'corrected'} = 1;
}

sub is_corrected {
    my ($self) = @_;
    return $self->{'corrected'};
}

sub add_seq_to_DS {
    my ($self, %args) = @_;
    if ( exists $args{'seqid'} && exists $args{'cluster'} ) {
        $self->{'DS_seqs'}->{ $args{'seqid'} } = $args{'cluster'};
    }
}

sub add_seq_to_SD {
    my ($self, %args) = @_;
    if ( exists $args{'seqid'} && exists $args{'cluster'} ) {
        $self->{'SD_seqs'}->{ $args{'seqid'} } = $args{'cluster'};
    }
}

sub SD_seqs {
    my ($self) = @_;
    return $self->{'SD_seqs'};
}

sub DS_seqs {
    my ($self) = @_;
    return $self->{'DS_seqs'};
}

# get cluster which seq was DS in
sub DS_seq_to_cluster {
    my ($self, $value) = @_;
    if ( exists $self->{'DS_seqs'}->{$value} ) {
        return $self->{'DS_seqs'}->{$value};
    }
}

# get cluster which seq was SD in
sub SD_seq_to_cluster {
    my ($self, $value) = @_;
    if ( exists $self->{'SD_seqs'}->{$value} ) {
        return $self->{'SD_seqs'}->{$value};
    }
}

# a sequence is misplaced only if it occurs in both the SD and DS lists
# so we only care about those in both. This allows us to track which
# seqs and which clusters are wrong
sub tabulate_SD_and_DS {
    my ($self) = @_;
    
    foreach my $seq (keys %{ $self->SD_seqs } ) {
        # my $cluster = $self->DS_seq_to_cluster($seq);
        # if (defined $cluster) {
        if (exists $self->{'DS_seqs'}->{$seq}) {
            print join("\t", 'MOVED', $seq, $self->SD_seq_to_cluster($seq),
             $self->DS_seq_to_cluster($seq) ), "\n";
        }
    }
}

1;
