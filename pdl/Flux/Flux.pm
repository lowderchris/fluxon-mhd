=head1 NAME

Flux - MHD simulation with fluxons

=head1 SYNOPSIS

use PDL;
use Flux;

$world = new Flux::World;

...

=head1 DESCRIPTION

FLUX is the Field Line Universal relaXer; it is a resistivity-free
magnetohydrodynamic simulator using the quasi-lagrangian fluxon
approach.  Fluxons are discretized magnetic field lines; they interact
at a distance via the magnetic curvature and pressure forces.

At the moment, FLUX is "merely" a topology-conserving force-free field
solver.  As development proceeds, MHS, QSMHD, and ultimately full MHD
will be incorporated into the code.

The Flux module is a perl interface to the underlying 'C' engine.  It
defines several objects that are part of the FLUX millieu:

=over 3

=item Flux::World

An entire simulation arena.  Each Flux::World object contains "global"
variables specific to that simulation.

=item Flux::Concentration

A boundary condition anchor

=item Flux::Fluxon

A discretized field line, implemented (in C) as a linked list of vertices,
connected by piecewise-linear flux elements ("fluxels").

=item Flux::Vertex

A single element of a fluxon.

=back

=head1 AUTHOR

FLUX is copyright (c) 2004 by Craig DeForest. Development was funded
under grants from NASA's Living With a Star program.  The code is
distributable under the terms of the Gnu Public License ("GPL"),
version 2.  You should have received a copy of the GPL with FLUX, in the
file named "COPYING".  

FLUX comes with NO WARRANTY of any kind.

Comments, questions, gripes, or kudos should be directed to
"deforest@boulder.swri.edu".

=cut

use Flux::World;
use Flux::Fluxon;
use Flux::Vertex;
use Flux::Concentration;


1;
