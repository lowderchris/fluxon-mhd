=head1 NAME

Flux::Vertex - fluxon vertex / fluxel

=head1 SYNOPSIS

  use PDL;
  use Flux;
  $world = new Flux::World;
  
  $fluxon = $world->fluxon(19);
  $vertex = $fluxon->vertex(4);
  $vertices = $fluxon->vertices;

=head1 Methods

=cut

BEGIN {
package Flux::Vertex;

require Exporter;
require DynaLoader;
@ISA = qw( Exporter DynaLoader );
@EXPORT = qw( );

bootstrap Flux::Vertex;
}

package Flux::Vertex;

use overload '""' => \&stringify;

=pod

=head2 stringify

=for ref 

Produce a string summary of a vertex ("" overload)

=cut

sub stringify {
    &_stringify(shift);
}



=pod

=head2 prev

=for ref

Return the previous vertex on the same fluxon

=head2 next

=for ref

Return the next vertex on the same fluxon

=cut

# Implemented in Vertex.xs



=pod

=head2 neighbors

=for ref

Return the neighbors of this vertex, as a perl list of VERTEX objects.

=head2 nearby

=for ref

Return the inverse-neighbor list of this vertex, as a perl list of VERTEX objects.

=cut

sub neighbors {
  my $this = shift;
  return @{_adjacent($this,0)};     # in Vertex.xs
}

sub nearby {
  my $this = shift; 
  return @{_adjacent($this,1)};     # in Vertex.xs
}

=pod

=head2 hull

=for ref

Return the hull of a vertex, in 2-D.

Calculates the neighborhood and 2-D projected hull, and returns the 
coordinates as a PDL.  The output is a 5xN PDL.  The columns are:
[2d-x, 2d-y, hull-x, hull-y, open]
where 2d-x and 2d-y are the projected coordinates of the neighbor points,
and hull-x and hull-y are the coordinates of the hull vertices.
"open" is a flag indicating whether the associated vertex is "open" (ie
nonexistent)

=cut

# Implemented in Vertex.xs

=head2 proj_neighbors

=for usage

    $xyl = $vertex->proj_neighbors(0); # just neighbor vertices
    $xyl = $vertex->proj_neighbors(1); # all vertices in the World

=for ref

Return just the neighbors, projected into the perpendicular plane of the 
vertex's segment using the fancy projection.  The output is a 3xN PDL.
The columns are: [2d-x,2d-y,label] where 'label' is the label of the other vertex.

Takes a single parameter - a global flag to indicate that all vertices, not just the
precalculated neighbors, are to be included in the calculation.

=cut

# Implemented in Vertex.xs

1;
