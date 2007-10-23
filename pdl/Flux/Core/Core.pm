=head1 NAME

Flux::Core - cheesy dynamic linker

=head1 SYNOPSIS

  use Flux::Core;

Flux::Core only holds a single global variable, $Flux::Core::FLUX,
which is itself a pointer to the FLUX "C" symbol table struct.  That's
a kludge to achieve cross-module linking of the same library.  Loading
Flux::Core links the FLUX C language library into the current running
PDL, and generates and populates a global C struct that contains
function pointers to all the entry points.  Then other modules can
access the library via the C struct, which is stored in Perl's address
space (and hence accessible to everyone).  Otherwise, each module
would link its own copy of the FLUX libraries, to the detriment of
all.

All the action happens in the BOOT: section of Core.xs.

=head1 Methods

-none-

=cut

  package Flux::Core;
  require DynaLoader;
  @ISA = qw/DynaLoader/;
  bootstrap Flux::Core;


1;
