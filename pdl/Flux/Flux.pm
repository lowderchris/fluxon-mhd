
use Flux::World;
use Flux::Fluxon;
use Flux::Vertex;
use Flux::Concentration;

BEGIN {
    package Flux;
    require DynaLoader;
    @ISA = qw(DynaLoader);
    bootstrap Flux;
}


=head1 NAME

Flux - the Field Line Universal relaXer (MHD without resistivity); v. 1.0

=head1 SYNOPSIS

  # (from within perl)
  use PDL;
  use Flux;
  $world = new Flux::World;
  <...>

=head1 DESCRIPTION

Flux is the name of the perl module that is used to interface to FLUX,
the Field Line Universal relaXer; it is a resistivity-free
magnetohydrodynamic simulator using the quasi-lagrangian fluxon
approach.  Fluxons are discretized magnetic field lines; they interact
at a distance via the magnetic curvature and pressure forces.

This is the man page for version 1.0 of Flux, the first publically
released version.  At the moment, FLUX is "merely" a
topology-conserving force-free field solver.  As development proceeds,
MHS, QSMHD, and ultimately full MHD will be incorporated into the
code.

The Flux module is a perl interface to the underlying 'C' engine.  It
defines several objects that are part of the FLUX millieu:

=over 3

=item Flux::World

An entire simulation arena.  Each Flux::World object contains "global"
variables specific to that simulation, and has access methods to 
get flux concentrations and fluxons in the simulation.

=item Flux::Concentration

A boundary condition anchor - each fluxon must begin and end at a flux
concentration.  Flux concentrations have locations, flux values, and 
types.  In practice flux concentrations from toy models typically have 
one fluxon apiece, though in empirically driven models they are anticipated
to require more.

=item Flux::Fluxon

A discretized field line, implemented (in C) as a linked list of vertices,
connected by piecewise-linear flux elements ("fluxels").  Each fluxon has a 
source and a sink flux concentration (Flux::Concentration), and is composed 
of a list of vertices (Flux::Vertex) that contain the actual spatial information.

=item Flux::Vertex

A single element of a fluxon.  A vertex contains location and connectivity information
to the next and previous vertices along the same fluxon.

=back

Each of those object types has its own perl module and its own man page (e.g. say
C<man Flux::World> to find out about what you can do with world objects).

VERSION

This is Flux version 1.1 - part of the FLUX 1.1 release.

=head1 GETTING STARTED

If you are reading this as a man page, you are most of the way there.
FLUX is a Perl module that relies on the Perl Data Language (PDL) to
work, as well as on the fluxlib.a library (written in C).  If you are
not acquainted with PDL you will want to make sure it is installed, by
entering the "pdl" command at your UNIX prompt (or if, God help you,
you are trying to do MHD simulations from within Microsoft Windows,
you should stop and install Fedora now).

The Flux module itself handles basic direct manipulation tasks; to run
a simulation, you will need some additional code that is included in the 
"pdl/PDL" directory of the main source code distribution.  Some of that code might
be included in this module in a later version.  The ".pdl" files in the PDL 
directory are intended to work with the PDL autoloader.  Your PDLLIB environment
variable, or the @PDLLIB global list within pdl, determines where the autoloader
searches (see L<PDL::AutoLoader> for details).

In particular, Flux does NOT contain an outermost loop for conducting 
relaxations. That is accomplished in perl space and is regarded as simple
enough that you can hand-roll your own as necessary.  There is a template in 
"pdl/PDL/simple_relaxer.pdl", which you can invoke with 
C<simple_relaxer($world, $int, $glob, $n, $opt)>, where C<$world> is a Flux::World
object, $int is an interactive flag (usually 0), $global performs rigorous neighbor
checking (at great expense in speed), and $n is the desired step number.  The $opt
parameter is an options hash; consult the source code for details.

Here is a sequence of commands that should help you get your first FLUX relaxation to work:
  
  use Flux;
  $fluxdir = "/usr/local/src/flux";       # or whatever directory you use
  push(@PDLLIB,"+$fluxdir/pdl/PDL");
  $a = read_world('/usr/local/src/flux/menagerie/potential.flux');
  $a->forces('f_pressure_equi','f_curvature','f_vertex');
  $a->b_flag(0);
  $dt = 0.3;
  simple_relaxer($a,0,0,1000);

The "pdl/menagerie" directory contains several other sample relaxation
description files that you can use to test your code.  See L<Flux::World> 
for details.

=head1 FILE FORMAT

Flux worlds are stored and transferred in an ASCII format that is parsed and generated
by the underlying "C" library.  Here is a description of the file format.

Flux world files are keyword/value data sets in standard ASCII form.  Lines that are
all whitespace or that begin with the '#' character are comments and are ignored.

Each world description file begins with some global information coupled with description of
the contents of the world or with a description of how the world or its boundary conditions are to
change.  The file format is intended to be useful both to store single-frame relaxation sets
and also to store and transfer time-dependent datasets.

Each line begins with a keyword followed by some parameters, separated
by whitespace.  Most of the characters of the keyword are ignored; these are marked with square brackets,
as are all other optional parts of the line.  Quantities to be replaced with numbers or strings are
surrounded with angle-brackets.  

All entities receive numeric labels that should be unique.  

=head2 Global description lines

Here are descriptions of line types that describe global changes to the simulation or 
delimit parts of the file.

=over 3

=item FR[AME] <Frame_no>

Marks the start of a new frame in time-dependent data

=item FI[NISHED] 

Marks the end of the data for a time step.

=item G[LOBAL] 

Global lines are used to set up the simulation, and there are several:

=over 3

=item GLOBAL F[ORCE] <forcename>,<forcename>,...

Sets which code force laws are to be used for the relaxation.  (Check the top of physics.c to 
find a list of the currently valid ones.  As of version 1.0, the B-normalized (DEPRECATED) forces
are correct, but the non-normalized forces are not yet validated.  Most relaxations require at
least three forces:  magnetic pressure (e.g. 'f_pressure_equi'), magnetic curvature (e.g. 'f_curvature'),
and a vertex pseudoforce (e.g. 'f_vertex').  If magnetic field measurements are required, then 
one of the 'b' forces must be included -- these are force laws that do not directly apply a force,
just calculate the magnetic field at each fluxel center.  

=item GLOBAL BO[UNDARY] <type> <par1> <par2> <par3> <par4> <par5> <par6>

Sets a global boundary condition.  Ultimately multiple global boundaries may be enabled but for now
at most one may be imposed.  The <type> code is one of "NONE", "PLANE", "SPHERE", or CYL".  For NONE, 
no additional parameters are required.  For PLANE, the first three numeric parameters are the origin of
the plane, and the second three are a vector normal to the plane.  For SPHERE, the first three numeric 
parameters are the origin of the sphere, and the fourth parameter is the radius; the other two are
ignored.  For CYLINDRICAL boundaries, only a single parameter -- the major radius of the cylinder -- is 
used.  The cylinder is always placed at the origin, with axis along the Z direction.

=item GLOBAL BF[LAG] <flag>

Sets whether B-normalized or direct forces are to be used.  If the flag is zero, then the forces
are B-normalized; if it is 1, then normal forces are in use.  (This affects the scaling of spatial
steps during the relaxation).

=item GLOBAL STATE <state>

Indicates the currennt state of the fluxon World as saved.  Not used in version 1.0.

=back

=back

The following types of line are used to describe the topology within the simulation arena:

=over 3

=item N[EW] <label> <x> <y> <z> <flux> [<type>]

Create a new flux concentration at the specified location with the specified amount of flux.  if
the type field is not specified, then the concentration is assumed to have fixed location.

=item D[ELETE] <label> 

Causes the labelled flux concentration to be deleted (mostly useful for time dependent data)

=item M[OVE] <label> <x> <y> <z>

Moves the corresponding flux concentration to the specified location.

=item L[INE] <fc_label1> <fc_label2> <fluxon_label> <flux>

Creates a fluxon connecting from the first flux concentration (which should have positive
flux) to the second (which should have negative flux). 

=item V[ERTEX] <FL_lab> <pos> <x> <y> <z>

Adds a vertex to the fluxon with label FL_lab, at the <pos> position.  Counting starts at 
zero, but assigning to position zero makes no sense as the 0th vertex is associated with the
originating flux concentration.  To assign to the end, use either a huge position number or
-2.  (0 is the originating location and 1 is the end location; you don't want to write to either
of them)


=cut


=back

=head1 METHODS

The Flux module itself is a shell that contains no external methods at
all -- see the individual object modules for details on available
methods.

There are several low-level internal methods that are used for
accessing FLUX data structures; you do not want to use them unless you
Really Know What You're Doing.

=cut

package Flux;

{
    our $typecodes = {
	links =>    1,
	vertex =>   2,
	fluxon =>   3,
        world  =>   4,
	concentration=>5
	};
    
    our $methods = {
	num    => [\&_rnum,  \&_wnum  ],
	long   => [\&_rlong, \&_wlong ],
	vector => [\&_rvec,  \&_wvec  ],
	vec    => [\&_rvec,  \&_wvec  ],
	Vertex => [\&_rvertex, undef  ],
	Fluxon => [\&_rfluxon, undef  ],
	FluxonList => [\&_rfluxonlist, undef],
	Concentration => [ \&_rconcentration, undef ],
        Neighbors => [sub{ _rdumblist("Flux::Vertex",@_) },
		      undef
		      ],
        Nearby => [sub{ _rdumblist("Flux::Vertex",@_) },
		      undef
		      ]
    };

    our $codes = { 
	vertex => {
	    line=>        [1, 'Line'],
	    prev=>        [2, 'Vertex'],
	    next=>        [3, 'Vertex'],
	    x=>           [4, 'vector'],
	    neighbors=>   [5, 'Neighbors'],
	    nearby=>      [6, 'Nearby'],
	    scr=>         [7, 'vector'],
	    r=>           [8, 'num'],
	    a=>           [9, 'num'],
	    b_vec=>      [10, 'vec'],
	    b_mag=>      [11, 'num'],
	    f_v=>        [12, 'vec'],
	    f_s=>        [13, 'vec'],
	    f_t=>        [14, 'vec'],
	    f_s_tot=>    [15, 'vec'],
	    f_v_tot=>    [16, 'num'],
	    r_v=>        [17, 'num'],
	    r_s=>        [18, 'num'],
	    r_cl=>       [19, 'num'],
	    label=>      [20, 'long'],
	    links_sum=>  [21, 'num'],
	    links_n=>    [22, 'long'],
	    links_up=>   [23, 'Vertex'],
	    links_left=> [24, 'Vertex'],
	    links_right=>[25, 'Vertex'],
	    links=>      [26, undef],
	    energy=>     [27, 'num']
	    },
	fluxon => {
	    flux=>	          [1,'num'],
	    label=>               [2,'long'],
	    start=>               [3,'Vertex'],
	    end=>                 [4,'Vertex'],
	    v_ct=>                [5,'long'],
	    all_links=>           [6,undef],
	    all_links_sum=>       [7,'num'],
	    all_links_n=>         [8,'long'],
	    all_links_up=>        [9,'Fluxon'],
	    all_links_left=>      [10,'Fluxon'],
	    all_links_right=>     [11,'Fluxon'],
	    start_links=>         [12,undef],
	    start_links_sum=>     [13,'num'],
	    start_links_n=>       [14,'long'],
	    start_links_up=>      [15,'Fluxon'],
	    start_links_left=>    [16,'Fluxon'],
	    start_links_right=>   [17,'Fluxon'],
	    end_links=>           [18,undef],
	    end_links_sum=>       [19,'num'],
	    end_links_n=>         [20,'long'],
	    end_links_up=>        [21,'Fluxon'],
	    end_links_left=>      [22,'Fluxon'],
	    end_links_right=>     [23,'Fluxon'],
	    fc_start =>           [24,'Concentration'],
	    fc_end =>             [25,'Concentration'],
	    },
	world => {
	    frame_number=>       [1,'long'],
	    state =>              [2,'long'],
	    concentrations =>    [3,'Concentration'],
	    lines =>             [4,'Fluxon'],
	    vertices =>          [5,'Vertex'],
	    photosphere =>       [6,undef],
	    image =>             [7,'Vertex'],
	    image2 =>            [8,'Vertex'],
	    locale_radius =>     [9,'num'],
	    fc_ob =>            [10,undef],
	    fc_oe =>            [11,undef],
	    fc_pb =>            [12,undef],
	    fc_pe =>            [13,undef],
	    verbosity =>        [14,'long'],
	    forces =>           [15,undef],
	    scale_b_power=>     [16,'num'],
	    scale_d_power=>     [17,'num'],
	    scale_s_power=>     [18,'num'],
	    scale_ds_power=>    [19,'num'],
	    refct=>             [20,'long'],
	    rc_funcs=>          [21,undef]
	    },
	concentration => {
	    world=>		 [1,'World'],
	    flux=>               [2,'num'],
	    label=>              [3,'long'],
	    lines=>              [4,'Fluxon'],
	    links=>	 	 [5,undef],
	    links_sum=>          [6,'num'],
	    links_n=>            [7,'long'],
	    links_up=>           [8,'Concentration'],
	    links_left=>         [9,'Concentration'],
	    links_right=>        [10,'Concentration'],
	    x=>                  [11,'vector'],
	    locale_radius=>      [12,'num'],
	}
    };



    ##############################
    # Assemble an ordering for next iteration through tied hashes -- sort the links
    # in each tied hash and generate a hash that links by name.
    # Kludgey but it works.

    our $ordering = {};

    foreach my $type(keys %$codes) {
	$ordering->{$type} = {};
	my $thash = $codes->{$type};

	my @order = (
		     sort
		        { $thash->{$a}->[0] <=> $thash->{$b}->[0] } 
		        keys %$thash
		     );

	for my $i(0..$#order) {
	    $ordering->{$type}->{$order[$i]} = $order[$i+1];
	}
    }

}    
	

=pod

=head2 r_val - read a value from a structure

=for usage

   FLUX::r_val( structcode, valcode, typestr );

=for ref

The structcode is the numeric code associated with the structure type; it should be one of the
values from the global hash ref $Flux::typecodes.  The valcode is a numeric type corresponding
to structure field, and typestr is the type string in the global $Flux::codes hash.  If it is a
string, it is used to look up the reader in the global $Flux::methods hash, or if it is an array
ref containing two code refs then they are executed directly.

=cut

sub r_val {
    my( $me, $struct, $field, $type ) = @_;
    my $reader;

    if(ref $type) {
	die "r_val: Can't handle ref types yet\n";
    } else {
	$reader = $Flux::methods->{$type}->[0];
    }

    unless( defined($reader)) {
	print "r_val: no read method for type $type ($struct,$field) - returning undef\n";
	return undef;
    }

    &{$reader}( $me, $struct, $field );
}

=pod

=head2 w_val - write a value to a structure

=for usage

    FLUX::w_val( $structcode, $valcode, $typestr, $value );

=for ref

The structcode is the numeric code associated with the structure type; is should be one of the
values from the global hash ref $Flux::typecodes.  The valcode is a numeric type corresponding
to the structure field, and typestr is the type string in the global $Flux::codes hash.  If it is a 
string, it is used to look up the reader in the global $Flux::methods hash, or if it is an array
ref containing two code refs then they are executed directly.  The value is the value to put in the 
structure field.

=cut

sub w_val {
    my( $me, $struct, $field, $type, $val ) = @_;
    my $writer;

    if(ref $type) {
	die "w_val: can't handle ref types yet\n";
    } else {
	$writer = $Flux::methods->{$type}->[1];
    }

    die "w_val: write method is undefined for type $type ...\n"
	unless defined($writer);

    &{$writer}( $me, $struct, $field, $val );
}

=pod

=head2 tree - traverse a tree structure and return all its leaves in a list ref

=for usage

    @list = FLUX::tree($root, $lname);

=for ref

    The lname is the prefix for the link fields in the structure (e.g. "all_links" for the fluxon
tree).  The routine could/should be subclassed as appropriate to pass in individual links fields.

=cut

sub tree {
    my $me = shift;
    my $lname = shift;
    my $whole_tree = shift;

    while($whole_tree && defined($me->{$lname."_up"})) {
	$me = $me->{$lname."_up"};
    }
    
    unless($me->{$lname."_n"}>1) {
	return $me;
    }

    my @out = ();
    my $l = $me->{$lname."_left"};
    my $r = $me->{$lname."_right"};
    
    push(@out, tree($l,$lname)) if(defined $l);
    push(@out, $me);
    push(@out, tree($r,$lname)) if(defined $r);
    @out;
}

=pod

=head1 AUTHOR

FLUX is copyright (c) 2004-2005 by Craig DeForest. Development was funded
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

$Flux::file_versions = Flux::file_versions();

1;
