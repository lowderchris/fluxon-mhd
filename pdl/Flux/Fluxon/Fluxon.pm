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
use overload '""' => \&_stringify;



1;
