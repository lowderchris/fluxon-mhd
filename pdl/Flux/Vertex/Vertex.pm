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


1;
