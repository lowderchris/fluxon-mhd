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



1;
