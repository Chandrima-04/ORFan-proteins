package Cluster;

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use ClusterUtils;

sub new {
    my ( $class, %args ) = @_;

    my $object = {
        'id'      => q{},
        'size'    => 0,
        'members' => {},
        'known'   => {},
    };
    bless( $object, $class );

    if    ( defined $args{'id'} )  { $object->{'id'} = $args{'id'}; }
    elsif ( defined $args{'-id'} ) { $object->{'id'} = $args{'-id'}; }
    else { die "must specify an ID for the Cluster!\n"; }

    return $object;
}

sub id {
    my ( $self, $id ) = @_;
    if ($id) { $self->{'id'} = $id; }
    else     { return $self->{'id'}; }
}

sub size {
    my ( $self, $size ) = @_;
    if ($size) { $self->{'size'} = $size; }
    else       { return $self->{'size'}; }
}

sub add_member {
    my ( $self, @args ) = @_;

    foreach my $one (@args) {
        if ( !( ref($one) eq 'Bio::Seq' ) ) {
            die "attempting to add something that's not a Bio::Seq\n";
        }
        else {
            $self->{'members'}->{ $one->display_id() } = $one;
        }
    }

    $self->{'size'} = scalar keys %{ $self->{'members'} };
    if ( $self->size() == 0 ) {
        die( "failed to add a seq to ", $self->id(), "\n" );
    }
    return $self;
}

sub get_member {
    my ( $self, @args ) = @_;

    my $seqobj = $self->{'members'}->{ $args[0] };

    if   ($seqobj) { return $seqobj; }
    else           { return; }
}

sub get_all_members {
    my ($self) = @_;

    return $self->{'members'};
}

sub add_known {
    my ( $self, @args ) = @_;

    if ( @args != 2 ) {
        die "can't pass more than 2 arguments to add_known\n";
    }
    my ( $id, $hit ) = @args;

    if ( $self->is_known($id) ) {
        warn "$id is already listed as known (hit $hit)\n";
    }
    $self->{'known'}->{$id} = $hit;

    return $self;
}

sub is_known {
    my ( $self, @args ) = @_;

    if (@args) {
        my $known_count = 0;
        foreach my $id (@args) {
            $known_count++ if defined( $self->{'known'}->{$id} );
        }
        return $known_count;
    }
    else {
        my $boolean = 0;
        if ( scalar %{ $self->{'known'} } > 0 ) { $boolean = 1; }
        return $boolean;
    }
}


# args:
# -fh       : give the filehandle for a fasta file or stream
# -normalize: boolean. remove lcl| from beginning of display IDs
sub add_members_from_fh {
    my ( $self, %args ) = @_;

    if ( !defined $args{'-fh'} ) {
        warn "ERROR: no filehandle. Can't add sequences!\n";
        return;
    }

    # grab seqs from the seqs filehandle and add them to the cluster
    my $in = Bio::SeqIO->new(
        -fh     => $args{'-fh'},
        -format => 'fasta',
    );

    while ( my $seq = $in->next_seq() ) {

        # get rid of lcl| prefix in seqID
        if ($args{'-normalize'}) {
            normalize_seqid($seq);
        }
        
        $self->add_member($seq);
        if ( $self->size() == 0 ) {
            die( $seq->display_id(), " didn't get added\n" );
        }
    }

    # check for an empty cluster before adding the new cluster
    if ( $self->size() == 0 ) {
        die( "You're trying to add an empty cluster: ", 
             $self->id(), "\n" );
    }

    return $self;
}


1;
