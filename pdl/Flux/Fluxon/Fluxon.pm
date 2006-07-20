=head1 NAME

Flux::Fluxon - discretized field line object

=head1 SYNOPSIS

  use PDL;
  use Flux;

  $world = new Flux::World;
  $fluxon = $world->fluxon(19);
  print $fluxon;

=head1 DESCRIPTION

Flux::Fluxon objects are the perl representation of fluxons within the
MHD simulator.  They are represented as tied hashes, to keep the interface
as flexible as possible.

VERSION

This is Fluxon.pm version 1.1, part of the FLUX 1.1 release.

=head1 METHODS

=cut

BEGIN {
package Flux::Fluxon;

require Exporter;
require DynaLoader;
@ISA = qw( Exporter DynaLoader Flux);
@EXPORT = qw(  ) ;

bootstrap Flux::Fluxon;
}

package Flux::Fluxon;
use overload '""' => \&stringify;


=pod

=head2 new_from_ptr 

=cut

sub new_from_ptr {
    my $class = shift;
    my $ptr = shift;
    my %hash;
    tie %hash,"Flux::Fluxon",$ptr;
    my $me= bless(\%hash,$class);
    _inc_world_refct($me);
    $me;
}


=pod

=head2 stringify

=for ref

Generate a string summary of a fluxon; overloaded with "".

=cut

sub stringify {
    &_stringify(shift);
}


=pod

=head2 vertex

=for usage

  $vertex = $fluxon->vertex(2);

=for ref

Retrieve the nth vertex from a particular fluxon (n begins at 0)

=cut

# Implemented in Fluxon.xs



=pod

=head2 vertices

=for usage
  
  @vertices = $fluxon->vertices();

=for ref

Retrieve all the vertices in a fluxon, as a perl list.

=cut

sub vertices { 
  my $me = shift;
  my $ct = $me->{v_ct};
  return map { vertex($me,$_) } 0..$ct-1;
}


=pod

=head2 polyline

=for usage

  $polyline = $fluxon->polyline;

=for ref

Return all the vertex locations in a fluxon, as a 3xn PDL.

(mainly useful for visualizations).  the 0th dim is, of course, (x,y,z).

=cut

# Implemented in Fluxon.xs

=pod

=head2 bfield

=for usage
  
  $bfield = $fluxon->bfield;

=for ref

Return the B vector at each vertex location in the simulation.

=cut

# Implemented in Fluxon.xs

=pod

=head2 dump_vecs

=for usage

  $stuff = $fluxon->dump_vecs;

=for ref

Returns a xN PDL containing, in each row:

=over 3

=item cols  0- 2: vertex location

=item cols  3- 5: B field vector

=item cols  6- 8: following segment partial force

=item cols  9-11: vertex partial force
    
=item col   12:   sum-of-magnitudes for the segment force components

=item col   13:   sum-of-magnitudes for the vertex force components

=item col   14:   r_s  - projected neighborhood radius for segment forces

=item col   15:   r_v  - projected neighborhood radius for vertex forces

=item col   16:   r_cl - closest neighbor approach projected radius

=back

=cut

# Implemented in Fluxon.xs

sub DESTROY {

  # DESTROY gets called twice -- once for the tied hash and once for the underlying object
  # (which is a scalar ref, not a hash ref). Ignore the destruction for the tied hash.

  eval ' my $a = ${$_[0]}; $a; ';
  return if($@);
  
  _dec_refct_destroy_world( $_[0] );

}

######################################################################
# TIED INTERFACE
# Mostly relies on the general utility functions in Flux....

sub TIEHASH {
    my $class = shift;
    my $ptr = shift;
    my $me = \$ptr;
    bless($me,$class);
#    _inc_world_refct($me);  
    return $me;
}

sub FETCH {
    my($me, $field)=@_;
    my $code = $Flux::codes->{fluxon}->{$field};

    return undef unless defined($code);
    Flux::r_val( $me, $Flux::typecodes->{fluxon}, @$code[0..1] );
}

sub STORE {
    my($me, $field,$val) = @_;
    my $code = $Flux::codes->{fluxon}->{$field};
    return undef unless defined($code);
    Flux::w_val( $me, $Flux::typecodes->{fluxon}, @$code[0..1], $val );
}

sub DELETE {
    print STDERR "Warning: can't delete fields from a tied FLUXON hash\n";
    return undef;
}

sub CLEAR {
    print STDERR "Warning: can't clear a tied FLUXON hash\n";
    return undef;
}

sub FIRSTKEY {
    return "flux";
}

sub NEXTKEY {
    my ($class,$prev) = @_;
    return $Flux::ordering->{fluxon}->{$prev};
    
}

sub SCALAR {
    _stringify(@_);
}


1;
