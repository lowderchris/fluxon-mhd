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

=head1 METHODS

=cut

BEGIN {
package Flux::Fluxon;

require Exporter;
require DynaLoader;
@ISA = qw( Exporter DynaLoader );
@EXPORT = qw(  ) ;

bootstrap Flux::Fluxon;
}

package Flux::Fluxon;
use overload '""' => \&stringify;


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
    @{_vertices(shift)};   # meat is in Fluxon.xs; unwrap the array ref.
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

=back

=cut

# Implemented in Fluxon.xs

1;
