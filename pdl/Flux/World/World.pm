=head1 NAME

Flux::World - Fluxon MHD simulation arena

=head1 SYNOPSIS

  use PDL;
  use Flux;

  $world = new Flux::World;

  $world = Flux::World->str2world($string);
  $world = Flux::World->read_world($filename);

  $str = $world->string();
  $ok  = $world->write_world($filename);

  $world->step($dt);
  $lines = $world->vectors;

=head1 DESCRIPTION

Flux objects are the perl interface to FLUX MHD models.  A Flux::World
object is the interface to a global aspects of the simulation.
A Flux::World object is just a pointer off into hyperspace.  


=head2 AUTHOR

Copyright 2004-2008 Craig DeForest.  You may distribute this file 
under the terms of the Gnu Public License (GPL), version  2. 
You should have received a copy of the GPL with this file.  If not,
you may retrieve it from "http://www.gnu.org".

=head2 VERSION

This file is part of the FLUX 2.2 release (22-Nov-2008)

=head1 FUNCTIONS

=cut

BEGIN {

use PDL::Graphics::Gnuplot 2.006_001;

package Flux::World;
use PDL;

require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader Flux);
@EXPORT = qw( read_world write_world str2world );

bootstrap Flux::World;
}

package Flux::World;
use overload '""' => \&_stringify;
use strict;
use warnings;

=pod

=head2 new

=for usage

$world = new Flux::World;

=cut

sub new {
    return _new_world();
}

=head2 read_world

=for usage

$world = read_world($file);

=for ref

Read in a fluxon structure from a file

You supply the file name; the fluxon gets read into memory.
The file is a "standard" FLUX ASCII fluxon description file.

If the file ends with a '.gz' extension, it is automatically
gunzipped on-the-fly.

=cut

sub read_world {
  shift if(substr($_[0],0,length(__PACKAGE__)) eq __PACKAGE__);
  my $file = $_[0];
  unless($file =~ m/.gz$/) {
      return _read_world(@_);
  } else {
      my $tmpfile = "/tmp/world-$$-tmp.flux";
      `gunzip < $file > $tmpfile`;
      my $w = _read_world($tmpfile);
      unlink $tmpfile;
      return $w;
  }
}


=pod

=head2 str2world

=for usage

$world = str2world($str)
$world = str2world(@lines)
$world = str2world([@lines])

=for ref

Convert string or collection of lines to a FLUX World

This allows you to convert a perl string to a fluxon simulation
structure.  It is implemented cheesily: the 'c' code uses streams to
interpret the data, so your strings are written to a temporary file
that is then snarfed back in.

=cut

sub str2world {
  shift if(substr($_[0],0,length(__PACKAGE__)) eq __PACKAGE__);

    if(@_ == 1 && ref $_[0] eq 'ARRAY') {
	@_ = @{$_[0]};
    }
    my $fname = "/tmp/fluxtmp-$$";
    open FLUXSYSTEM_TMP,">$fname" || die "Can't open tmp file $fname\n";
    print FLUXSYSTEM_TMP join("\n",@_,"");
    close FLUXSYSTEM_TMP;
    my $world = read_world($fname);
    unlink $fname;
    $world;
}



=pod

=head2 write_world

=for usage

    write_world($world,$file);
    $world->write_world($file);

=for ref

Write a fluxon structure out to a file

Writes an ASCII description file to the filename you give.

If the file name ends in ".gz", the file is automatically gzipped on-the-fly.

=cut

sub write_world {
    my $world = shift;
    my $fname = shift;
    
    unless($fname =~ m/\.gz$/) {
	_write_world($world,$fname);
    } else {
	my $tmpfname = "/tmp/fluxtmp-$$";
	eval {unlink $tmpfname;};
	my $wwo = _write_world($world, $tmpfname);
	`gzip $tmpfname`;
	`mv ${tmpfname}.gz $fname`;
	return $wwo;
    }
}
    

=pod

=head2 string

=for ref

Convert a fluxon model to an ASCII string or (in list context) list.

Works mighty cheesily, by using the C code to write the string out to 
a temporary file, then snarfing the temp file back in.


=cut

sub string {
    my $me = shift;
    
    my $tmpfile = "/tmp/fluxtmp-$$";
    eval {unlink $tmpfile;};
    write_world($me,$tmpfile);
    open FLUXSYSTEM_TMPFILE,"<$tmpfile" || die "Couldn't open temp. file\n";
    my @lines = <FLUXSYSTEM_TMPFILE>;

    return join("",@lines) unless wantarray;
    return @lines;
}


=pod

=head2 _stringify

=for ref

Produce a short summary string describing a fluxon model.  This routine is used to convert
a World to an ASCII string whenever it is used in string context.


=cut
    
sub _stringify {
    my $me = shift;

    my @s = $me->string;
    my @lines = grep(m/^\s*LINE/,@s);
    my @vertices = grep(m/^\s*VERTEX/,@s);
    my @bounds = ( sub {"None. "}
		   ,sub {"Plane; origin: $_[0],$_[1],$_[2]; normal: $_[3],$_[4],$_[5]";}
		   ,sub {"Sphere; origin: $_[0],$_[1],$_[2]; radius: $_[3]"}
		   ,sub {"Cylinder origin: $_[0],$_[1],$_[2], normal: $_[3],$_[4],$_[5], radius:".sqrt($_[3]**2+$_[4]**2+$_[5]**2)}
		   );
    my @ph = $me->photosphere;
    my @ph2 = $me->photosphere2;
    my $btype = $ph[6];
    my $btype2 = $ph2[6];

    my $scalehash = $me->step_scales();


    no warnings;    # shut up about undefined elements while printing
    return "Fluxon geometry object; ".(@lines+0)." fluxons, ".(@vertices+0)." vertices\n"
	."\tGlobal boundary: ".&{$bounds[$btype]}(@ph)."\n"
	."\tGlobal b2: ".&{$bounds[$btype2]}(@ph2)."\n"
	."\tForce scaling powers: ".join("; ",(map { "$_=$scalehash->{$_}" } keys %$scalehash))."\n\n"
	;
    use warnings;

}

=pod

=head2 summary

=for ref

Produce a longer summary string describing a fluxon model.  This works by dumping the 
value of each hash element on the tied-hash side.

=cut

sub summary {
    my $me = shift;
    my $s = "";
    for my $k(keys %$me) {
	$s .= sprintf( "%15s %s\n", $k.":", $me->{$k} );
    }
    $s;
}

=pod step_scales

=for usage

    $href = $world->step_scales;  # pulls step scales out of the world
    $world->step_scales($href);   # puts them in

=for ref

(Deprecated; use the hash interface instead, via the scale_FOO_power fields)

Starting with FLUX 1.1, there is fine grained control over the power laws used for
calculating the size of the steps taken by individual vertices.  The force laws all
calculate force per unit length along the fluxon, so length compensation is needed; 
also, the largest steps should in general be taken where the field is weak (and the
forces are weak) so some scaling is required.  The step law is:

    dx = F * dt * d^pd * B^pb * S^ps * DS^pds

where d is the harmonic mean of relevant neighbor distances, B is the
magnetic field strength, S is the stiffness coefficient of the force
calculation (vector sum / magnitude sum  of the component forces), and 
DS is the difference-stiffness relating the point's force to its neighbors' forces.

These global values are communicated to Perl via hash ref that is passed into and
out of the WORLD. 

=cut

# Implemented in World.xs



=pod

=head2 fix_proximity

=for usage

    $baddies = $world->fix_proximity($thresh);

=for ref

Adds additional vertices as needed to ensure that no fluxels get too close.
The $thresh is the ratio of fluxel length / closest approach.  Fluxels that
violate that criterion get divided by the addition of more vertices.
You get back a count of the number of fluxels that required subdivision.

You have to have primed the world with at least one timestep or 
update_neighbors call before you do this -- otherwise horrible things might
happen.

=cut

# Implemented in World.xs




=pod

=head2 fix_curvature

=for usage

    $num_changed = $world->fix_curvature($thresh_high, $thresh_low);
    $num_changed = $world->fix_curvature($thresh_high);
    
=for ref

Adds and subtracts vertices as needed to ensure that no vertices have too
much or too little curvature.  The $thresh_high is the maximum bend angle, in
radians, that is allowed at each vertex.  $thresh_low is the minimum
bend angle in radians, though no vertex is removed unless it meets certain
proximity criteria as well. You get back a count of the number of fluxels
different from the previous number.

You have to have primed the world with at least one timestep or
update_neighbors call before you do this -- otherwise horrible things might
happen.

=cut

sub fix_curvature { 
    my $me = shift ;
    my $thresh_high = shift ;
    my $thresh_low = shift ;
    _fix_curvature ($me, $thresh_high || 0, $thresh_low || 0 );
}

=pod

=head2 update_neighbors

=for usage

    $world->update_neighbors($globalflag);

=for ref

Update the neighbor trees between adjacent fluxels.  If $globalflag is
nonzero, then all previous neighbor information is discarded, causing
a global neighbor search for each fluxel.  The global algorithm is
O(n^2).  If $globalflag is zero, then a next-nearest-neighbors
algorithm is used, which is O(n).  

The first time you call this, you had better set $globalflag=1, or
horrible things will happen.  After that, you will want to set $globalflag=0
so that you don't die of ennui before the simulation completes.

=cut

# Implemented in World.xs


=pod

=head2 verbosity

=for usage

    $v = $world->verbosity();
    $world->verbosity($v);


=for ref

(Deprecated - use the 'verbosity' field in the tied hash instead)

Controls the level of, well, verbosity associated with operations on $world.
High values send insane amounts of stuff to the console.

Some approximate meanings for each verbosity level:

=over 3

=item 0 

Be as quiet as possible

=item 1

Describe global items and iffy conditions

=item 2

Describe activity on a per-fluxon level

=item 3 

Describe activity on a per-vertex/per-fluxel level

=item 4

Describe activity within each vertex operation -- notably neighbor searches.

=back

=cut

# Implemented in World.xs


=pod

=head2 fluxon_ids

=for usage

  @list = $world->fluxon_ids;
  $fluxon_count = $world->fluxon_ids;

=for ref

Return a list of all fluxon labels in the simulation

Walks through the World all-fluxon tree and gives you the label numbers 
of all the fluxons, as a perl list.

=cut

sub fluxon_ids {
  my $ref = _fluxon_ids(@_); # Call _fluxons in World.xs
  @$ref;
}


=pod

=head2 concentrations

=for usage

    @concentrations = $world->concentrations;

=for ref

Returns a list of all flux concentrations in the world.

=cut

sub _conc_helper {
    my $me = shift;
    my $depth = shift||0;
    my @a;

    if(exists($me->{links_left})) {
	push(@a, _conc_helper($me->{links_left},$depth+1));
    }

    push(@a, $me);

    if(exists($me->{links_right})) {
	push(@a, _conc_helper($me->{links_right},$depth+1));
    }

    return @a;
}

   
sub concentrations {
    my $me = shift;
    my @list =_conc_helper($me->{concentrations}); 
    print "_conc_helper returned ".(0+@list)." elements....\n" if $me->{verbosity};
    return @list;
}

=pod

=head2 concentration

=for usage

 $conc = $world->concentration($cid);

=for ref

Return the flux concentration wtih the given ID, as a Flux::Concentration object.

=cut

# Implemented in World.xs

=pod

=head2 fluxon

=for usage
  
  $fluxon = $world->fluxon($fid);

=for ref

Return the fluxon with the given ID, as a Flux::Fluxon object.  

=cut

# fluxon is implemented in World.xs

=pod

=head2 fluxons

=for usage
 
  @fluxons = $world->fluxons(@fids);

=for ref

Return the fluxon(s) with the given IDs, as a list of Flux::Fluxon objects.  If you
give no arguments, you get all of them.

=cut

sub fluxons {
    my $world = shift;
    if(!@_) {

      push(@_,$world->fluxon_ids);
    }
    my $id;
    my @fluxons;
    while(defined ($id=shift)) {
	push(@fluxons, &fluxon($world,$id)) unless($id<0 && $id>=-11);
    }
    
    return @fluxons;
}

=pod

=head2 vertex

=for usage

    $vertex = $world->vertex($vid);

=for ref

Return the vertex with the given ID, as a Flux::Vertex object.

=cut

# vertex is implemented in World.xs    


=head2 vertex_ids

=for usage

  @list = $world->vertex_ids;
  $vertex_count = $world->vertex_ids;

=for ref

Return a list of all vertex labels in the simulation

=cut

sub vertex_ids {
  my $ref = _vertex_ids(@_);
  @$ref;
}



=pod

=head2 vertices

=for usage
 
 @vertices = $world->vertices(@v_ids)

=for ref

Return the vertices with the given IDs, as a list of Flux::Vertex objects.  if you
give no arguments, you get all of them (sorted in numerical order).

=cut

sub vertices {
    my $world = shift;
    my @a = @_;
    if(!@a){
	push(@a,@{$world->_vertex_ids});
    }
    my $id;
    my @vertices;
    while(defined ($id=shift @a)) {
	push(@vertices, $world->vertex($id));
    }
    return @vertices;
}

=pod

=head2 new_concentration - generate a new flux concentration

=for usage

    $fc = $world->new_concentration( $where, $flux, $label, $end );

=for ref

Adds a new flux concentration to the world at the given location --
$where must be a 3-PDL or something convertible into one.  $flux is the
amount of magnetic flux to be stored in the flux concentration.  If $label
is nonzero then it will be used as the label field of the new concentration, 
otherwise a unique label will be autochosen.  If $end is present then it
should be an end-condition string; otherwise the worldwide default end condition
is used.

=cut

sub new_concentration {
    my $w = shift;
    my $x = pdl(shift);
    my $flux = shift;
    my $label = shift;
    my $endsv = shift;
    return _new_concentration($w,$x,$flux,$label,$endsv); # _new_concentration in World.xs
}


=pod

=head2 emerge - generate two flux concentrations connected by a fluxon

=for usage

    $source = $world->emerge($source_loc, $sink_loc, $flux, $vertices);
    $source = $world->emerge($source_loc, $sink_loc, $flux, $n);

=for ref

This is a convenience routine for producing two flux concentrations and a single-fluxon 
loop at a desired location.  You specify the locations of the two concentrations (as 3-PDLs)
and the flux (currently ignored of course but included as a placeholder for later).  If
you hand in a 3xN PDL it is treated as a vertex location list.  If you hand in a scalar 
(or 1 PDL) it is treated as a number of nontrivial vertices to include, and the fluxon is
shaped as close to a straight line as possible consistent with the photosphere in the world.

The flux concentrations are line-tied with the default boundary.

=cut

sub emerge {
    my($world, $src, $snk, $flux, $vertices) = @_;

    barf "emerge needs a world"        unless(ref $world eq 'Flux::World');
    barf "emerge needs a 3-PDL source" unless(ref $src eq 'PDL' and $src->ndims==1 and $src->dim(0)==3);
    barf "emerge needs a 3-PDL sink"   unless(ref $snk eq 'PDL' and $snk->ndims==1 and $snk->dim(0)==3);
    
    $flux = 1 unless($flux);
    
    barf "vertices must be a scalar or 3xN" if(ref $vertices eq 'PDL' and 
					       $vertices->dim(0) != 3
					       );

    my $fc0 = new_concentration($world, $src, 1, 0, undef);
    print "emerge: after first new_concentration, ref count is $world->{refct}\n" if($world->{verbosity});
    my $fc1 = new_concentration($world, $snk, -1, 0, undef);
    print "emerge: after 2nd new_concentration, ref count is $world->{refct}\n" if($world->{verbosity});

    my $fl;

    if( ref $vertices eq 'PDL' && $vertices->dim(0)==3) {

	$fl = $fc0->new_fluxon($fc1, 1, 0, $vertices);
	print "emerge: after fluxon origination, ref count is $world->{refct}\n" if($world->{verbosity});

    } else {
	$fl = $fc0->new_fluxon($fc1, 1, 0);
	print "emerge: after fluxon origination, ref count is $world->{refct} (nonspecified vertex case)\n" if($world->{verbosity});

	if("$vertices") {
	    my $nv = ( pdl(0)+$vertices )->at(0);

	    print "emerge: generating vertices ($vertices)\n" if($world->{verbosity});

	    my $svec = ($snk - $src)/($nv+1);
	    my $sep = (($snk - $src) * ($snk - $src))->sumover->sqrt;
	    my $ph = $world->{photosphere};
	    my $vertex = $fl->{start};
	    for my $i(1..$nv){
		my $loc = $svec * $i + $src;

		if($ph->{type}==1) {
		    ### Planar photosphere -- if it's below the photosphere, move it up to the surface
		    ### plus 10% of the separation between the footpoints.

		    my $vertcoord = (($loc - $ph->{origin}) * ($ph->{normal}))->sumover;
		    if($vertcoord < 0) {
			$loc -= $vertcoord * $ph->{normal};
			$loc += $sep * 0.1 * $ph->{normal};
		    }

		} elsif($ph->{type}==2) {
		    ### Spherical photosphere -- if it's below the photosphere, move it outward to 
		    ### the surface plus 1% of the separation between the footpoints.
		    
		    my $radius = (($loc - $ph->{origin}) * ($loc - $ph->{origin}))->sumover->sqrt;
		    if($radius < $ph->{normal}->((0))) {
			my $rhat = ($loc - $ph->{origin}) / $radius;
			my $disp = $ph->{normal}->((0)) + 0.1 * $sep - $radius;
			$loc += $disp * $rhat;
		    }
		}
		$vertex = $vertex->add_vertex_after( $loc );
	    }
	}
    }
    return $fc0;
}
    

=pod

=head2 forces - read or set the force laws used for a world

=for usage
 
  print $world->forces;
  $world->forces('t_curvature','t_vertex')

=for ref

(Deprecated; use the "forces" field in the tied hash instead)

Return the names of the forces that are queued up to act on each fluxel
in the simulation, or set them.  The forces have predefined names
that are present at the top of physics.c.  Unknown forces are referred to
by hex value, but that way lies madness unless you are the linker (and
I don't mean you, John).

Like all the original access methods, "forces" is deprecated -- you ought
to be using the more regularized hash interface instead (the relevant
field is $a->{forces}).  

=cut

sub forces {
    my $me = shift;
    if(@_==0) {
	return @{_forces($me)};    # Retrieve forces (in World.xs)
    } else {
	my $i = 0;
	my $s;
	while($s=shift) {
	    _set_force($me,$i++,$s); # Set forces one at a time (in World.xs)
	}
	_set_force($me,$i,'');
    }
}

=pod

=head2 reconnection - read or set the reconnection criteria

=for usage

  print $world->reconnection;
  $world->reconnection('rc_a_ad2',[3.14159/4, 10 * 3.14159/4], ...);

=for ref

(Deprecated; use the "rc_funcs" field in the tied hash instead)

The reconnection criteria accept numeric parameters in an array ref.
Check the documentation array or the source code to find out what they
mean.  Setting the reconnection parameters doesn't make reconnection
happen -- for that, you must call "reconnect".  If you add more than one
criterion, just string 'em together.

=cut

sub reconnection {
  my $me = shift;
  if(@_==0) {
    return @{_rcfuncs($me)};
  } else {
    my $i=0;
    while(@_) {
      my $name = shift;
      my $params = shift;
      die "error in _set_rc\n" unless(_set_rc($me, $i++, $name, $params));
    }
    _set_rc($me,$i++,"",[]);
  }
}


=pod

=head2 boundary

=for usage

 $phot = pdl($world->photosphere);
 $world->boundary([$phot->list],$type);

=for ref

(Deprecated; use the photosphere element of the tied hash instead)

With no additional arguments, returns a seven-element perl list.  Elements
0-6 are the numeric values of the boundary parameters, and element 7 is
the type code for the boundary.  If no boundary is in use, then the empty
list is returned.

If you pass in a list ref, then the contents should be a 
seven-element list containing the 
(origin, normal-vector) coordinates of the photospheric plane, followed
by a numeric type code, or a zero-element list indicating that there is no 
photospheric plane.  If you
pass in a list ref, it is used to set the photospheric plane parameters.
The list ref is copied into the boundary parameter array internally.  The 
$type code should be the numeric type code for the boundary condition you want:
0 is none, 1 is planar, 2 is spherical, and 3 is cylindrical.  More might be
added later.  This is something of a kludge at the moment -- string parsing should
be in here (as in the force law parameters).

=cut

sub boundary { 
  return @{Flux::_photosphere(@_)};
}
sub photosphere {
    return @{Flux::_photosphere(@_)};
}
sub photosphere2 {
    return @{Flux::_photosphere($_[0],$_[1]||undef,$_[2]||undef,32)};
}


=pod

=head2 update_force

=for usage

$world->update_force;
$world->update_force($globalflag);

=for ref

Walks through the world and updates the forces on each vertex.  You should 
update the neighbors before the first time you run update_force.
(the C side of this is called "world_update_mag" for historical reasons).
The $globalflag determines whether the neighbor determination should be 
global (nonzero; slow) or local (zero; default; much faster).

=cut

# Implemented in World.xs

=pod

=head2 relax_step

=for usage

$world->relax_step($dt);

=for ref

Updates positions of the vertices after the forces have been
calculated.  Takes one step toward relaxation.  Warning -- if $dt is
too big you'll be in trouble, guv!  $dt is scaled so that 1.0 should
carry each vertex all the way to relaxation if it were the only vertex
in motion.  To avoid various sorts of awfulness, you want to keep it
to 0.3-0.5 in the early stages of most relaxations.

=cut

# Implemented in World.xs


=pod

=head2 fw_stats

=for usage

print $a->fw_stats->{n_av};

=for ref

Produces a statistical overview of the simulation, in a hash:  average force
per unit length on each vertex, average unsigned-total force on each vertex,
and average number of neighbors per vertex.  The answer comes back in a perl hash ref.
Mnemonic: Flux World stats.

=cut

# Implemented in World.xs
*fw_stats = \&_stats;

=pod

=head2 render

=for usage

 $world->render($interactive, $range, $opt)

=for ref

Produces a simple 3-D plot of all the field lines in a World, using PDL::Graphics::Gnuplot.

The C<$interactive> argument is a flag indicating whether the event loop should
be run to position and angle the output.  C<$range>, if present, is a 3x2 PDL containing
the (minimum, maximum) corners of the 3-cube to render.  C<$opt> is an options hash.
Currently useful options are:

=over 3

=item dev (default 'wxt')

Gnuplot device type to use for plotting.  Needs to be a terminal that accepts the 'dashed' term option.  Suggested terminals that should be common across different operating systems are (static files:) 'pngcairo', 'pdfcairo', 'postscript', 'svg', (interactive:) 'wxt', 'x11'.

=item rgb

If present, specifies that all lines should have this color (3-PDL)

=item prgb

If present, specifies that all points should have this color (3-PDL)

=item dip_detector

If this is set to a 3-PDL, then dips in the magnetic field are
rendered in this color regardless of what color they would otherwise have.
Dips are segments adjacent to a vertex that is lower than both its neighbors.

=item dip_detector2

If this is set to a 3-PDL, then dips in the magnetic field are
rendered in this color regardless of what color they would otherwise have.
In this case, the entire dip is recorded out to the altitude of the lower
"lip" of the dip, not just the segments adjacent to the bottom of the dip.

=item no_dip

If this is set, it should be a hash ref whose keys are the labels of
particular fluxons that should *not* be dip-detected by the
dip_detector or dip_detector2 options.  This is useful for masking out
a particular locus of the sim that should be dip-detected.  

=item rgb_fluxons

If present, this should be a hash ref whose keys are fluxon id numbers or the string
'default', and whose values are either:

=over 3

=item a 3-PDL (RGB triplet) for the fluxon.  

This sets the whole fluxon to the given color

=item a pair of 3-PDLs with the colors for the north and south poles respectively.  The colors will be interpolated.

The fluxon color is interpolated between these two values.

=item a triplet of 3=PDLs with the colors for north, center, and south poles respectively.  The colors will be interpolated.

=back

=item RGB_CODE

If present, contains a code ref that calculates the RGB values of
individual line segments in the local context of render.  (Not for the
faint of heart!)

=item PRGB_CODE

If present, contains a code ref in the same style as RGB_CODE that
calculates the RGB values of individual points in the local context of
render. (Not for the faint of heart!)

=item concentrations

Flag indicating whether or not to draw flux concentrations separately at the
end of field lines.  They are rendered as large points. (Default is 1).

=item rgb_fcs

If present, this should be a hash ref whose keys are flux concentration labels
and whose values are 3-PDLs with the color in which the concentration is to be 
rendered.  (By default, source concentrations are blue, sink concentrations are
red, unknown concentrations are white, and errors are green).

=item points

Flag indicating whether or not to draw vertices separately from field lines,
using a Points object.  (Default is 1).  

=item psize

Flag indicating the size, in screen pixels, of the points as rendered.
(Default is 1; 0 is not allowed, but fractional values are)

=item linewidth

Flag indicating the width of the LineStrip objects used to render the fluxons.
(Default is 1; 2-3 is useful for generating figures)

=item label

Flag indicating whether to indicate the numeric label of every vertex.
(Default is 0).  If label is an array ref, then it is treated as a
list of fluxons to label -- all other fluxons are NOT labeled.

=item label_fluxons

Flag indicating whether to label individual fluxons (at both endpoints).

=item neighbors

Flag indicating whether to render the closest approach of each fluxel's neighbors 
in the perpendicular plane to the fluxel. (default is 0)

=item hull

Flag indicating whether to render the hull of each fluxel in the perpendicular plane of the
fluxel.  (Default is 0)

=item nscale

The scaling factor of the trace-lines and neighbor points (default is
0.25).  1.0 produces nice intersection diagrams showing the
approximate segmentation of space into the neighborhoods of each
fluxon, but is quite busy.  The default is a compromise between loss
of detail and separation of neighborhood tracings.

=back

=cut

use PDL;
use PDL::Graphics::Gnuplot;
use PDL::NiceSlice;
use PDL::Options;

*render_lines = \&render;
our $window;

sub render {
    my $w = shift;
    my $opt=shift // {};
    my $dev = $opt->{dev} // 'wxt';
    my $gpwin = shift // $opt->{window} // $window // ($window=gpwin($dev,size=>[9,9],dashed=>0));

    $gpwin->options(trid=>1);

    my (@rgb,@prgb);
    print "Defining RGB..." if($Flux::debug);
    my @fluxons = ();
    for my $id($w->fluxon_ids) {
	push(@fluxons,$id) if($id>0 || $id<-11);
    }

    my $fid;
    my @poly = map { print "$_..." if($Flux::debug); $w->fluxon($_)->polyline } @fluxons;


    if($opt->{'RGB_CODE'}) {
      ### RGB_CODE - just execute the code string
      eval $opt->{'RGB_CODE'};

    } elsif($opt->{'rgb_fluxons'}) {
      ### RGB_FLUXONS - specify color per-fluxon
      @rgb = ();

      print "rgb_fluxons..." if($Flux::debug);

      for my $i(0..$#fluxons) {

	my $spec = $opt->{'rgb_fluxons'}->{$fluxons[$i]+0};
	unless (ref $spec eq 'PDL' || (ref $spec eq 'ARRAY' and @$spec > 0)) {
	  $spec = $opt->{'rgb_fluxons'}->{default};
	  unless (ref $spec eq 'PDL' || (ref $spec eq 'ARRAY' and @$spec > 0)) {
	    $spec = pdl(1,1,1);
	  }
	}

	# Convert 3 x n PDLs to lists of 3-PDLs
	if(ref $spec eq 'PDL' and $spec->dims == 2 ) {
	  $spec = [ dog $spec ];
	}

	print "fluxon $i: spec is $spec...\n" if($Flux::debug);

	if(ref $spec eq 'PDL') {

	  push(@rgb, $spec * ones($poly[$i]));

	} elsif(ref $spec eq 'ARRAY') {

	  if( @$spec == 1 ) {
	    push(@rgb, $spec->[0] * ones($poly[$i]));

	  } elsif( @$spec == 2 ) {
	    my $alpha = double yvals($poly[$i]) / (($poly[$i]->dim(1) - 1)||1);
	    my $beta = 1.0 - $alpha;
	    push(@rgb,  ($alpha * $spec->[0] +
			 $beta  * $spec->[1]
			 ));

	  } elsif( @$spec >= 3 ) {
	    my $alpha = double yvals($poly[$i]) / (($poly[$i]->dim(1) - 1)||1);
	    my $beta = 1.0-$alpha;
	    my $gamma = sin($alpha * 3.14159);
	    push(@rgb, ($alpha * $spec->[0] +
			$beta  * $spec->[2] +
			$gamma * $spec->[1]
			)
		 );
	  }

	} else {

	  ## Error -- should not happen ##
	  push(@rgb, pdl(1,0.25,0) * ones($poly[$i]));

	}
      }

    } elsif(defined $opt->{'rgb'}) {
      ### RGB - specify color for all fluxons at once
      @rgb = map { (ones($_) - (yvals($_)==0)) * $opt->{'rgb'} } @poly;

    } else {
      ### Default case - blue->red color scheme for all fluxons
      @rgb = map { 
	my $alpha = double yvals($_);
	$alpha /= max($alpha);
	    my $beta = 1.0 - $alpha;
	    my $gamma = sin($alpha*3.14159)**2;
	    my $prgb = $beta * pdl(1,0,0) + $alpha * pdl(0,0,1) + $gamma * pdl(0,1,0);
	    $prgb;
	} @poly;
    }

    print "Defining PRGB..." if($Flux::debug);
    if($opt->{'PRGB_CODE'}) {

      eval $opt->{'PRGB_CODE'};

    } elsif(defined $opt->{'prgb'}) {

      @prgb = map { (ones($_) - (yvals($_)==0)) * $opt->{'rgb'} } @poly;

    } else {

      @prgb = map { 
	my $alpha = double yvals($_);
	$alpha /= max($alpha);
	my $beta = 1.0 - $alpha;
	my $gamma = sin($alpha*3.14159)**2;
	my $prgb = $beta * pdl(1,0,0) + $alpha * pdl(0,0,1) + $gamma * pdl(0,1,0);
	$prgb;
      } @poly;
    }


    ##############################
    # Having defined the colors of everything, check for sophisticated dip detection and execute if necessary
    if(ref $opt->{'dip_detector2'} eq 'PDL') {
      print "Detecting dips...\n" if($Flux::debug);
      # Walk through the fluxons and identify any midpoint vertices that are lower than their 
      # neighbors on the same fluxon
      for my $i(0..$#poly) {
	next if((ref $opt->{'no_dip'} eq 'HASH') && ($opt->{'no_dip'}->{$fluxons[$i]}));
	my $poly = $poly[$i];
	my $stack = $poly->range([[2,-1],[2,0],[2,1]],[0,$poly->dim(1)],'e');
	my $mask =  ($stack->((1)) <= $stack->((0))) & ($stack->((1)) <= $stack->((2)));
	$mask->(0) .= 0;
	$mask->(-1) .= 0;

	# Now spread outward in the forward direction till we reach a local maximum
	my @dips = list which($mask);
	my @mw = ();

	for my $dip(@dips) {
	    my ($before, $after) = ($dip, $dip);
	    while($before >= 2 && $poly->(2,$before-1)>$poly->(2,$before)) {
		$before--;
	    }
	    while($after < $poly->dim(1) - 1 && $poly->(2,$after+1)>$poly->(2,$after) &&
		  ($after==$dip || $poly->(2,$after+1)<= $poly->(2,$before))) {
		$after++;
	    }
	    while($before < $dip-1 && $poly->(2,$before) > $poly->(2,$after)) {
		$before++;
	    }

	    push(@mw,$before..$after);
	}
	
	if(@mw) {
	    my $mw = pdl(@mw);
	    my $prgb = $prgb[$i];
	    if($mw->nelem) {
		my $prgb = $prgb[$i];
		$prgb->(:,$mw) .= $opt->{'dip_detector2'};
		my $rgb= $rgb[$i];
		$rgb->(:,$mw-1) .= $opt->{'dip_detector2'};
		$rgb->(:,$mw) .= $opt->{'dip_detector2'};
	    }
	}	
      }
    }


    ##############################
    # Having defined the colors of everything, check for simplistic dip detection and execute if necessary
    if(ref $opt->{'dip_detector'} eq 'PDL') {
      print "Detecting dips...\n" if($Flux::debug);
      # Walk through the fluxons and identify any midpoint vertices that are lower than their 
      # neighbors on the same fluxon
      for my $i(0..$#poly) {
	my $poly = $poly[$i];
	my $stack = $poly->range([[2,-1],[2,0],[2,1]],[0,$poly->dim(1)],'e');
	my $mask =  ($stack->((1)) <= $stack->((0))) & ($stack->((1)) <= $stack->((2)));
	$mask->(0) .= 0;
	$mask->(-1) .= 0;
	print "mask is $mask\n" if($Flux::debug);
	my $mw = which($mask);
	my $prgb = $prgb[$i];
	if($mw->nelem) {
	  $prgb->(:,$mw) .= $opt->{'dip_detector'};
	  my $rgb= $rgb[$i];
	  $rgb->(:,$mw) .= $opt->{'dip_detector'};
	  $rgb->(:,$mw+1) .= $opt->{'dip_detector'};
	}
	
      }
    }

    print "Defining polygons...\n" if($Flux::debug);
    my $poly = null->glue(1,@poly);
    print "a...\n" if($Flux::debug);
    my $rgb = null->glue(1,@rgb);
    print "b...\n" if($Flux::debug);
    my $prgb = null->glue(1,@prgb);

    $gpwin->options($opt->{popt}) if(ref($opt->{popt}) eq 'HASH');

    my @plot;

    if($opt->{points} || !defined($opt->{points})) {
	push @plot,{with=>'points',lc=>'rgb variable',pointsize=>($opt->{psize}||1)},$poly->using(0,1,2),$prgb->(-1:0)->mult(255,0)->shiftleft(pdl(16,8,0),0)->sumover;
    }

    ##############################
    # Generate line strips for each fluxon...
    my $rgbdex = 0;
    my @id = @fluxons;

    my $olabel = pdl($opt->{label}) if ( ref $opt->{label} eq 'ARRAY');

    print "Defining line strips...\n" if($Flux::debug);
    for my $i(0..$#id) {
	my $fp = $w->fluxon($id[$i])->polyline;

	push @plot,{with=>'lines',lc=>'rgb variable',lw=>($opt->{linewidth}||1)},$fp->using(0,1,2),$rgb[$i]->(-1:0)->mult(255,0)->shiftleft(pdl(16,8,0),0)->sumover;
	#DL3D could do this and the 'points' in one step with a switch between 'lines' or 'linespoints'? Yes, except for the prgb/rgb duplication

	##############################
	# Label fluxon endpoints if desired
	if($opt->{label_fluxons}) {
	    $gpwin->options(label=>["N-$id[$i]-N",at=>join(',',$fp->(:,(0) )->list),'front']);
	    $gpwin->options(label=>["S-$id[$i]-S",at=>join(',',$fp->(:,(-1))->list),'front']);
	}

	##############################
	# Label every vertex if desired
	if($opt->{label}) {

	    next if(ref $opt->{label} eq 'ARRAY' && 
		    which( $id[$i] == $olabel )->nelem != 1
		    );

		
	    my @labels = map { $_->id } $w->fluxon($id[$i])->vertices ;

	    if(ref $opt->{label} eq 'ARRAY') {
		
	    }
		    
	    label:for my $j(0..$#labels) {
		$gpwin->options(label=>[qq/"$labels[$j]"/,at=>join(',',$fp->(:,($j))->list),'front']);
	    }
	}


    }

    ##############################
    # Generate flux concentrations if desired

    if($opt->{'concentrations'} || !defined($opt->{'concentrations'})) {
	my @fc_points = ();
	my @fc_rgb = ();
	for my $fc( $w->concentrations ) {
	    if($fc->{label} < -99 || $fc->{label} > 0) {
		push( @fc_points, $fc->{x} );
		if($opt->{'rgb_fcs'} && defined $opt->{'rgb_fcs'}->{$fc->{label}}) {
		    push(@fc_rgb, $opt->{'rgb_fcs'}->{$fc->{label}});
		} else {
		    if($fc->{lines}) {
			if( $fc->{lines}->{fc_start}->{label} == $fc->{label} ) {
			    push(@fc_rgb, pdl(0,0.25,1));
			} elsif( $fc->{lines}->{fc_end}->{label} == $fc->{label} ) {
			    push(@fc_rgb, pdl(1,0.25,0));
			} else {
			    push(@fc_rgb, pdl(0.25,1,0.25));
			}
		    } else {
			push(@fc_rgb, pdl(1,1,1));
		    }
		}
	    }
	}
	if(@fc_points) {
	    my $fc_points = cat(@fc_points);
	    my $fc_rgb = cat(@fc_rgb);
	    push @plot,{with=>'points',lc=>'rgb variable',pointsize=>2},$fc_points->using(0,1,2),$fc_rgb->mult(255,0)->shiftleft(pdl(16,8,0),0)->sumover;
	}
    }
	    

    print "Neighbors?\n" if($Flux::debug);
    my $nscale = $opt->{'nscale'} || 0.25;

    if($opt->{'neighbors'}){
	my @neighbors;

	for my $v (map { $_->vertices } $w->fluxons) {
	    next unless($v->next && ($v->id<-9 || $v->id>0));
	    my $xcen = 0.5 * ($v->x + $v->next->x);

	    my $pm = $v->projmatrix;

	    my $neigh = $v->proj_neighbors; 

	    for my $i(0..$neigh->dim(1)-1) {
		my $x = zeroes(3);
		$x->(0:1) .= $neigh->(0:1,($i));
		$x *= $nscale;
		my $xx = Flux::World::_vec_mmult_3d($pm->list,$x->list);
		$xx += $xcen;
		push(@neighbors,$xx);
	    }
	}
	
	my $fp = cat(@neighbors);
	push @plot,{with=>'points',pointsize=>2},$fp->using(0,1,2);

    }

    if($opt->{'hull'}) {

	my $hullrgb = defined($opt->{'hullrgb'}) ? $opt->{'hullrgb'} : pdl(0.3,0.3,0);
	print "hullrgb: $hullrgb...\n" if($Flux::debug);

	my $hullopen = defined($opt->{'hullopen'}) ? $opt->{'hullopen'} : 10;
	print "hullopen: $hullopen...\n" if($Flux::debug);

	my $zz = 0;
	for my $v( $w->vertices ) { ### map { $_->vertices } $w->fluxons) {
	    next unless($v->next && ($v->id<-9 || $v->id>0));

	    my $xcen = 0.5 * ($v->x + $v->next->x);
	    
	    my $pm = $v->projmatrix;
	    
	    my $hull = $v->hull(0);

	    my @hpoints = ();

	    for my $i(0..$hull->dim(1)-1) {

		my $rows = $hull->range([0,$i],[7,2],'p');
		if( ! ($rows->((4),:)->any) ) {
		    my $x = zeroes(3);
		    $x->(0:1) .= $rows->(2:3,(0));
		    $x *= $nscale;
		    my $xx = $x x $pm;
		    $xx += $xcen;
		    push(@hpoints,$xx);

		    $x->(0:1) .= $rows->(2:3,(1));
		    $x->(2) .= 0;
		    $x *= $nscale;
		    $xx = $x x $pm;
		    $xx += $xcen;
		    push(@hpoints,$xx);
		} elsif($rows->((4),0)) {
		    if(@hpoints) {
			my $fp = cat(@hpoints)->(:,(0),:);

			push @plot, {with=>'lines',lc=>[sprintf("#%06x",$hullrgb->mult(255,0)->shiftleft(pdl(16,8,0),0)->sumover)]},$fp->using(0,1,2);


			@hpoints = ();
		    }


		    my $x = zeroes(3);
		    $x->(0) .= $rows->((2),(1)) + $hullopen * cos($rows->((6),(0)));
		    $x->(1) .= $rows->((3),(1)) + $hullopen * sin($rows->((6),(0)));
		    $x *= $nscale;
		    my $xx = $x x $pm;
		    $xx += $xcen;
		    push(@hpoints,$xx);
		    
		    $x->(0:1) .= $rows->(2:3,(1));
		    $x->(2) .= 0;
		    $x *= $nscale;
		    $xx = $x x $pm;
		    $xx += $xcen;
		    push(@hpoints,$xx);
		} elsif($rows->((4),1)) {

		    my $x = zeroes(3);
		    $x->(0:1) .= $rows->(2:3,(0));
		    $x *= $nscale;
		    my $xx = $x x $pm;
		    $xx += $xcen;
		    push(@hpoints,$xx);

		    $x->(0) .= $rows->((2),(0)) + $hullopen * cos($rows->((5),(1)));
		    $x->(1) .= $rows->((3),(0)) + $hullopen * sin($rows->((5),(1)));
		    $x *= $nscale;
		    $xx = $x x $pm;
		    $xx += $xcen;
		    push(@hpoints,$xx);
		}
	    }

	    if(@hpoints) {

		my $fp = cat(@hpoints)->(:,(0),:);

		push @plot,{with=>'lines',lc=>[sprintf("#%06x",$hullrgb->mult(255,0)->shiftleft(pdl(16,8,0),0)->sumover)]},$fp->using(0,1,2);
	    }
	}    
    }

    @Flux::World::plotlist = @plot;

    unless($gpwin) {
	$gpwin = gpwin($dev,size=>[9,9]);
    }
    $gpwin->plot3d(@plot);
    $window = $gpwin;

}

=pod


=head2 polyline

=for ref

Dumps a giant polyline for all fluxons in the simulation.

=cut

sub polyline {
    my $w = shift;
    my $id_list = shift;
    my (@l,$l1);
    @l = map { $_->polyline } (defined($id_list) ? @$id_list : $w->fluxons);
    $l1 = $l[0];
    return $l1->glue(1,@l[1..$#l]);
}



=pod

=head2 bfield

Dumps the magnetic field value for all fluxons in the simulation

=cut

sub bfield {
    my $w = shift;
    my (@l,$l1);
    @l = map { $_->bfield } $w->fluxons;
    $l1 = $l[0];
    return $l1->glue(1,@l[1..$#l]);
}

=pod

=head2 energy

=for usage 

$w->energy

=for ref

Prints the global magnetic field energy (requires that you ran update_force
with one of the e_ forces enabled)

=cut

sub energy {

    my $w = shift;

    my $E = 0;
    
    for ($w->vertices) {
	$E += $_->{energy};
    }

    return $E;

}
	    

=pod

=head2 dump_vecs

=for ref

Dumps all vectors for all fluxons (see Flux::Fluxon::dump_vecs)

=cut

sub dump_vecs {
    my $w = shift;
    my (@l,$l1);
    @l = map { $_->dump_vecs } $w->fluxons;
    $l1 = $l[0];
    return $l1->glue(1,@l[1..$#l]);
}

=pod

=head2 closest_vertex

=for usage

$vertex = $world->closest_vertex(pdl($x,$y,$z), [$global, [$start_vertex]]);

=for ref

Return the closest vertex to a given spatial location.  At a minimum, you 
supply a 3-PDL containing the location.  The <code>$global</code> flag 
indicates whether a global search or a more direct linear search should
be performed.

The linear (non-global) search depends on the WORLD having a valid neighbor
calculation in place.

If a linear search is performed, it can be carried out on a previous vertex,
which will probably increase speed (especially if the previous vertex is close
to the desired location!).  If you feed in a <code>$start_vertex</code> and
the <code>$global</code> flag is 0, then the search starts from the given
start vertex.

You can call _closest_vertex instead and gain a little speed, but then you 
must include all arguments (even if they are simply undef).

=cut

sub closest_vertex {
    my $w = shift;
    my $x = shift;
    my $global = shift || 0;
    my $start_vertex = shift || undef;
    return _closest_vertex($w,$x,$global,$start_vertex);
}

=pod

=head2 closest_simplex

=for usage

@list = $world->closest_simplex(pdl($x,$y,$z), [$global, [$start_vertex]]);

=for ref

Return four vertices forming a tetrahedron that encloses the requested point.
If that is not possible (e.g. point is outside the outermost fluxon) then return
the nearest three.

The vertices come back as a list.  The parameters are the same as closest_vertex.

BEWARE: although the initial point is found by calling closest_vertex, subseqent 
stages of the simplex construction require neighbor initialization, so you can't just
set the global flag to 1 and expect everything to be OK on an unknown world -- it
must have appropriate neighbors already calculated.

=cut

sub closest_simplex {
    my $w = shift;
    my $x = shift;
    my $global = shift || 0;
    my $start_vertex = shift || undef;
    my $a = _closest_simplex($w,$x,$global,$start_vertex);
    if(ref $a) {
	return @$a;
    } else {
	return ();
    }
}

=pod

=head2 _hull_points

=for usage

    $hull = Flux::World::_hull_points($points);

=for ref

This is not a method at all, but rather a subroutine in the
Flux::World package, intended for testing and debugging the
minimum-hull routine.  It accepts a 2xN PDL containing points on the
plane, and returns a 6xM PDL containing hull information. On output,
the columns are as follows:

=over 3

=item 0,1

The coordinates of a neighbor point

=item 2,3

The coordinates of the corresponding leftward hull point if it exists, or the
(left,right) exit angles in radians if it doesn't.  

=item 4

A flag that is 0 if the hull point exists and 1 if it doesn't.

=item 5

The index of this neighbor in the original input array.

=back

Note that M <= N in all cases; rejected points are ignored.
The output is suitable to feeding into _plot_hull, below.

=cut

# implemented in World.xs

=pod

=head2 _plot_hull

=for usage

    $w = pgwin(dev=>"/xs", size=>[5,5]);
    Flux::World::_plot_hull($w, Flux::World::_hull_points($points) [, $points]);

=for ref

Interprets and plots a 6xN hull PDL returned by _hull_points.  If you
include the optional $points parameter, then those points are all plotted too.

If you want a particular title or scaling you should call $win->env beforehand and then 
hold $win.  $win is released on exit.

=cut

sub _plot_hull {
    my $win = shift;
    my $hull = shift;
    my $points = shift;

    # Plot all points in set, if presented
    if(defined $points) {
	$win->points($points->((0)),$points->((1)),{color=>2});
	$win->hold;
	for my $i(0..$points->dim(1)-1){
	    $win->text($i,$points->((0),($i)),$points->((1),($i)),{color=>2});
	}
    }

    # Plot all points in neighbors.
    $win->points($hull->((0)),$hull->((1)),{color=>3});
    $win->hold;
    for my $i(0..$hull->dim(1)-1){
	$win->text($hull->at(5,$i),$hull->at(0,$i),$hull->at(1,$i),{color=>3});
	$win->line(pdl(0,$hull->at(0,$i)),pdl(0,$hull->at(1,$i)),{color=>3});
    }

    print "Drawing hull....\n";
    # Draw all hull lines
    for my $i(0..$hull->dim(1)-1){
	my $r = $hull->range([2,$i],[3,2],'p');

	unless($r->at(2,0) || $r->at(2,1)) {
	    $win->line($r->((0)), $r->((1)),{color=>4});
	} elsif($r->at(2,0)) {
	    my $x1 = $r->at(0,1) + 2000 * cos($r->at(1,0));
	    my $y1 = $r->at(1,1) + 2000 * sin($r->at(1,0));
	    $win->line(pdl($r->at(0,1),$x1),pdl($r->at(1,1),$y1),{color=>5});
	} elsif($r->at(2,1)) {
	    my $x1 = $r->at(0,0) + 2000 * cos($r->at(0,1));
	    my $y1 = $r->at(1,0) + 2000 * sin($r->at(0,1));
	    $win->line(pdl($r->at(0,0),$x1),pdl($r->at(1,0),$y1),{color=>6});
	}
    }

    $win->points(0,0,{symsize=>3,charsize=>3,color=>1});

    $win->release;
}		


=head2 check

=for usage

    print $w->check;

=for ref

Performs some consistency checks on the World, and executes minor fixes where possible.  Returns 0 if the World
looked OK, 1 if it got fixed up, and -1 if it is damaged beyond simple repair.  

This is actually just a jump point into the C routine "world_check", in model.c.

=cut

## implemented in World.xs

=pod

=head2 interpolate_value - interpolate a numeric field from the VERTICES onto a specific spatial location

=for usage

    $val = $a->interpolate_value( $field_name, $loc, $global, $seed );
    $val = $a->interpolate_value( $field_name, $loc );

=for ref
    
Uses the initial string to identify a field type to interpolate, and returns an interpolated value for it. 
If the field type is a 'num', you get back the straight interpolation.  If it is a 'vec', you get back a 3-PDL 
containing all interpolated components.  Other types get undef.

If $global exists and is nonzero, then a global (rather than local) search is done.

=cut

sub interpolate_value {
    my $w = shift;
    my $field = shift;
    my $loc = shift;
    my $global = shift || 0;
    my $seed = shift || 0;
    my $offset;
    
    if((ref $w) ne 'Flux::World') {
	die "Flux::World::interpolate_value needs a Flux::World\n";
    }

    if( ((ref $loc) ne 'PDL') || ($loc->ndims != 1) || $loc->dim(0) != 3 ) {
	die "Flux::World::interpolate_value needs a 3-PDL location (no threading yet)\n";
    }
    
    if( $Flux::codes->{vertex}->{$field}->[1] =~ m/num/ ) {
	$offset = Flux::_ptr_offset( $Flux::typecodes->{vertex}, $Flux::codes->{vertex}->{$field}->[0] );
	return _interpolate_value($w, $loc, $global, $seed, $offset);
    }
    elsif( $Flux::codes->{vertex}->{$field}->[1] =~ m/vector/) {
	$offset = Flux::_ptr_offset( $Flux::typecodes->{vertex}, $Flux::codes->{vertex}->{$field}->[0] );
	return _interpolate_vector($w,$loc,$global,$seed,$offset);
    }	
    else {
	die "Flux::World::interpolate_value: field '$field' isn't a num or a vector...\n";
    }
}

## implemented in World.xs



sub DESTROY {
    Flux::destroy_sv( $_[0] );
}

######################################################################
# TIED INTERFACE
# Mostly relies on the general utility functions in Flux....

# TIEHASH not used much - the C side uses FLUX->sv_from_ptr to do the tying.
sub TIEHASH {
    my $class = shift;
    my $me = shift;
    return bless($me,$class);
}

sub EXISTS {
    my ($me, $field)=@_;
    return ($FLUX::codes->{world}->{$field});
}

sub FETCH {
    my($me, $field)=@_;
    my $code = $Flux::codes->{world}->{$field};
 
    return undef unless defined($code);
    
    Flux::r_val( $me, $Flux::typecodes->{world}, @$code[0..1] );
}

sub STORE {
    my($me, $field,$val) = @_;
    my $code = $Flux::codes->{world}->{$field};
    return undef unless defined($code);
    Flux::w_val( $me, $Flux::typecodes->{world}, @$code[0..1], $val );
}

sub DELETE {
    print STDERR "Warning: can't delete fields from a tied WORLD hash\n";
    return undef;
}

sub CLEAR {
    print STDERR "Warning: can't clear a tied WORLD hash\n";
    return undef;
}

sub FIRSTKEY {
    return "frame_number";
}

sub NEXTKEY {
    my ($class,$prev) = @_;
    return $Flux::ordering->{world}->{$prev};
    
}

sub SCALAR {
    _stringify(@_);
}


  
1;

