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

1;
