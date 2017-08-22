package ClusterGroup;

use strict;
use warnings;
use Storable qw(nstore retrieve);

# SYNOPSIS
# my $cluster_group = ClusterGroup->new( '-load' => $load );
# my @clusters = $cluster_group->get_all_clusters;


sub new {
    my ( $class, %args ) = @_;
    my $object;

    if ( $args{'-load'} ) {
        my $file = $args{'-load'};
        $object = retrieve($file);
        die "ERROR: unable to load from $file\n" unless defined $object;
    }
    else {
        $object = {
            'clusters' => {},    # key is cluster_id, value is Cluster obj
        };
    }
    bless( $object, $class );

    return $object;
}

sub add_cluster {
    my ( $self, $cluster ) = @_;

    my $id = $cluster->id;
    if ( exists $self->{'clusters'}->{$id} ) {
        warn "ERROR: cluster $id already exists!\n";
    }
    else {
        $self->{'clusters'}->{$id} = $cluster;
    }

    return $self;
}

sub get_cluster {
    my ( $self, %args ) = @_;

    my $id;
    if ( defined $args{'-id'} ) {
        $id = $args{'-id'};
    }

    if ( exists $self->{'clusters'}->{$id} ) {
        return $self->{'clusters'}->{$id};
    }
    else {
        warn "ERROR: no cluster with id $id!\n";
        return;
    }
}

sub get_all_cluster_ids {
    my ($self, %args) = @_;

    my @clusters;
    if (defined $args{'-sorted'}) {
        @clusters = sort { $a <=> $b } keys %{ $self->{'clusters'} };
    }
    else {
        @clusters = keys %{ $self->{'clusters'} };
    }

    return @clusters;    
}

# optional arg: -sorted
# $group->get_all_clusters(-sorted => 1);
# will return the clusters sorted numerically by clustername
sub get_all_clusters {
    my ($self, %args) = @_;

    my @clusters;
    if (defined $args{'-sorted'}) {
        my @cluster_names = sort { $a <=> $b } keys %{ $self->{'clusters'} };
        foreach my $name (@cluster_names) {
            push @clusters, $self->get_cluster('-id' => $name);
        }
    }
    else {
        @clusters = values %{ $self->{'clusters'} };
    }

    return @clusters;
}

sub count {
    my ($self) = @_;
    return scalar keys %{ $self->{'clusters'} };
}

sub save {
    my ( $self, $filename ) = @_;

    nstore( $self, $filename )
      or die "ERROR: couldn't save to file $filename: $!\n";

    return $self;
}

1;
