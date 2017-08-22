package Cluster;

use Bio::Seq;

sub new {
    my ($class, %args) = @_;
    
    my $object = {'id'      => q{},
		  'size'    => 0,
		  'members' => {},
		  'known'   => {},
                 };
    bless ($object, $class);
    
    if    (defined $args{'id'}  ) { $object->{'id'} = $args{'id'};  }
    elsif (defined $args{'-id'} ) { $object->{'id'} = $args{'-id'}; }
    else                          { die "must specify an ID for the Cluster!\n"; }

    
    return $object;
}

sub id { 
    my ($self, $id) = @_;
    if ($id) { $self->{'id'} = $id; }
    else { return $self->{'id'}; }
}

sub size { 
    my ($self, $size) = @_;
    if ($size) { $self->{'size'} = $size; }
    else { return $self->{'size'}; }
}

sub add_member {
    my ($self, @args) = @_;

    foreach my $one (@args) {
	if (!(ref($one) eq 'Bio::Seq')) {
	    die "attempting to add something that's not a Bio::Seq\n";
	}
	else {
	    $self->{'members'}->{$one->display_id()} = $one;
	}
    }

    $self->{'size'} = scalar keys %{$self->{'members'}};
    if ($self->size() == 0) {
	die ("failed to add a seq to ", $self->id(), "\n");
    }
    return $self;
}

sub get_member {
    my ($self, @args) = @_;

    my $seqobj = $self->{'members'}->{$args[0]};
    
    if ($seqobj) { return $seqobj; }
    else         { return; }
}

sub get_all_members {
    my ($self) = @_;
    
    my %members = %{ $self->{'members'} };

    return \%members;
}

sub add_known {
    my ($self, @args) = @_;

    if (@args != 2) {
	die "can't pass more than 2 arguments to add_known\n";
    }
    my ($id, $hit) = @args;

    if ($self->is_known($id)) {
	warn "$id is already listed as known (hit $hit)\n";
    }
    $self->{'known'}->{$id} = $hit;

    return $self;
}

sub is_known { 
    my ($self, @args) = @_;

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


1;
