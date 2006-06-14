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

VERSION

This is version 1.1 of World.pm - part of the 1.1 Flux release, 21-Oct-2005

=head1 FUNCTIONS

=cut

BEGIN {

use PDL::Graphics::TriD; 

package Flux::World;

require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader Flux);
@EXPORT = qw( read_world write_world str2world );

bootstrap Flux::World;
}

package Flux::World;
use overload '""' => \&_stringify;


sub new_from_ptr {
    my $class = shift;
    my $ptr = shift;
    my %hash;
    tie %hash,"Flux::World",$ptr;
    $hash{refct}=1;      # one reference (that we're returning)
    return bless(\%hash,$class);
}

=pod

=head2 read_world

=for usage

$world = read_world($file);

=for ref

Read in a fluxon structure from a file

You supply the file name; the fluxon gets read into memory.
The file is a "standard" FLUX ASCII fluxon description file.

=cut

sub read_world {
  shift if(substr($_[0],0,length(__PACKAGE__)) eq __PACKAGE__);
  _read_world(@_);
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
    @lines = <FLUXSYSTEM_TMPFILE>;

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
    my $s="";

    my @s = $me->string;
    my @lines = grep(m/^\s*LINE/,@s);
    my @vertices = grep(m/^\s*VERTEX/,@s);
    my @bounds = ( sub {"None. "}
		   ,sub {"Plane; origin: $_[0],$_[1],$_[2]; normal: $_[3],$_[4],$_[5]";}
		   ,sub {"Sphere; origin: $_[0],$_[1],$_[2]; radius: $_[3]"}
		   ,sub {"Cylinder (at origin, axis along Z); radius: $_[0]"}
		   );
    my @ph = $me->photosphere;
    $btype = $ph[6];

    my $scalehash = $me->step_scales();

    return "Fluxon geometry object; ".(@lines+0)." fluxons, ".(@vertices+0)." vertices\n"
	."\tGlobal boundary: ".&{$bounds[$btype]}(@ph)."\n"
	."\tForce scaling powers: ".join("; ",(map { "$_=$scalehash->{$_}" } keys %$scalehash))."\n\n"
	;

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

Adds additional vertices as needed to ensure that no vertices have too
much curvature.  The $thresh_high is the maximum bend angle, in
radians, that is allowed at each vertex.  $thresh_low is the minimum
bend angle in radians. You get back a count of the number of fluxels
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
      push(@fluxons, &fluxon($world,$id));
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

=head2 forces

=for usage
 
  print $world->forces;
  $world->forces('t_curvature','t_vertex')

=for ref

Return the names of the forces that are queued up to act on each fluxel
in the simulation, or set them.  The forces have predefined names
that are present at the top of physics.c.  Unknown forces are referred to
by hex value, but that way lies madness unless you are the linker (and
I don't mean you, John).

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

=head2 b_flag

=for usage
 
  print $world->b_flag;
  $world->b_flag(1);

=for ref

Returns the state of the B-normalization flag in the World.
When set to 0, the flag causes all the forces to be treated as if the 
local magnetic field magnitude were divided out.  This is appropriate
for the older forces (f_pressure_equi and f_curvature) but not
for the newer ones.  When set to 1, the flag causes all the forces
to be treated as, well, ordinary physical forces.

=cut

sub b_flag {
    my $me = shift;
    if(@_==0) {
	return $me->_b_flag;    # Retrieve forces (in World.xs)
    } else {
	$me->_set_b_flag(@_);
	return;
    }
}


=pod

=head2 boundary

=for usage

 $phot = pdl($world->photosphere);
 $world->boundary([$phot->list],$type);

=for ref

With no additional arguments, returns a seven-element perl list.  Elements
0-6 are the numeric values of the boundary parameters, and element 7 is
the type code for the boundary.  If no boundary is in use, then the empty
list is returned.

If you pass in a list ref, then the contents 

six-element list containing the 
(origin, normal-vector) coordinates of the photospheric plane, or a
zero-element list indicating that there is no photospheric plane.  If you
pass in a list ref, it is used to set the photospheric plane parameters.
The list ref is copied into the boundary parameter array internally.  The 
$type code should be the numeric type code for the boundary condition you want:
0 is none, 1 is planar, 2 is spherical, and 3 is cylindrical.  More might be
added later.  This is something of a kludge at the moment -- string parsing should
be in here (as in the force law parameters).

=cut

sub boundary { 
  return @{_photosphere(@_)};
}
sub photosphere {
    return @{_photosphere(@_)};
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

=head2 stats

=for usage

print $a->stats->{n_av};

=for ref

Produces a statistical overview of the simulation, in a hash:  average force
per unit length on each vertex, average unsigned-total force on each vertex,
and average number of neighbors per vertex.  The answer comes back in a perl hash ref.

=cut

# Implemented in World.xs
*stats = \&_stats;

=pod

=head2 render_lines

=for usage

 $world->render_lines($interactive, $range, $opt)

=for ref

Produces a simple 3-D plot of all the field lines in a World, using PDL::TriD.

The C<$interactive> argument is a flag indicating whether the eventt loop should
be run to position and angle the output.  C<$range>, if present, is a 3x2 PDL containing
the (minimum, maximum) corners of the 3-cube to render.  C<$opt> is an options hash.
Currently useful options are:

=over 3

=item rgb

If present, specifies that all lines should have this color (3-PDL)

=item prgb

If present, specifies that all points should have this color (3-PDL)

=item RGB_CODE

If present, contains a code ref that calculates the RGB values of
individual line segments in the local context of render.  (Not for the
faint of heart!)

=item PRGB_CODE

If present, contains a code ref in the same style as RGB_CODE that
calculates the RGB values of individual points in the local context of
render. (Not for the faint of heart!)

=item points

Flag indicating whether or not to draw vertices separately from field lines,
using a Points object.  (Default is 1).  

=item psize

Flag indicating the size, in screen pixels, of the points as rendered.
(Default is 4; 0 is not allowed)

=item linewidth

Flag indicating the width of the LineStrip objects used to render the fluxons.
(Default is 1; 2-3 is useful for generating figures)

=item label

Flag indicating whether to use a Label object to indicate the numeric label of every
vertex.  (Default is 0).

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
use PDL::Graphics::TriD;
use PDL::NiceSlice;

*render_lines = \&render;

sub render {
    my $w = shift;
    my $twiddle = shift;
    my $range=shift;
    my $opt=shift;

    if(!defined $opt && ref $range eq 'HASH') {
	$opt = $range;
	undef $range;
    }
    $opt= {} unless defined($opt);

    print "Releasing...\n" if($Flux::debug);
    release3d;
    print "foo...\n" if($Flux::debug);

    my $fid;
    my @poly = map { print "$_..." if($Flux::debug); $w->fluxon($_)->polyline } $w->fluxon_ids;

    my @rgb,@prgb;
    print "Defining RGB..." if($Flux::debug);
    if($opt->{'RGB_CODE'}) {
      eval $opt->{'RGB_CODE'};
    } elsif(defined $opt->{'rgb'}) {
      @rgb = map { (ones($_) - (yvals($_)==0)) * $opt->{'rgb'} } @poly;
    } else {
      @rgb = map { 
	my $alpha = double yvals($_);
	$alpha /= max($alpha);
	    my $beta = 1.0 - $alpha;
	    my $gamma = sin($alpha*3.14159);
	    my $prgb = $alpha * pdl(1,0,0) + $beta * pdl(0,0,1) + $gamma * pdl(0,1,0);
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
	    my $gamma = sin($alpha*3.14159);
	    my $prgb = $alpha * pdl(1,0,0) + $beta * pdl(0,0,1) + $gamma * pdl(0,1,0);
	    $prgb;
	} @poly;
    }

    print "Defining polygons...\n" if($Flux::debug);
    $poly = pdl(0,0,0)->glue(1,@poly);
    print "a...\n" if($Flux::debug);
    $rgb = pdl(0,0,0)->glue(1,@rgb);
    print "b...\n" if($Flux::debug);
    $prgb = pdl(0,0,0)->glue(1,@prgb);


    print "nokeeptwiddling3d...\n" if($Flux::debug);
    nokeeptwiddling3d;
    my ($boxmax,$boxmin);

    unless(defined($range)) {
	$range = cat($poly->mv(-1,0)->minimum, $poly->mv(-1,0)->maximum);
	my $rctr = $range->mv(-1,0)->average;
	my $rsize = ($range->(:,(1))-$range->(:,(0)))->max;
	$range->(:,(0)) .= $rctr - $rsize/2;
	$range->(:,(1)) .= $rctr + $rsize/2;
    }
	
    print "box definitions...\n" if($Flux::debug);
    $boxmax = $range->(:,(1));
    $boxmin = $range->(:,(0));

    my $shift =($boxmax-$boxmin)*0.05;
    $boxmin -= $shift;
    $boxmax += $shift;
    my $shift =($boxmax-$boxmin)*0.05;
    $boxmin -= $shift;
    $boxmax += $shift;

    print "Clipping polys...\n" if($Flux::debug);
    $pol2 = $poly->copy;
    $poly->((0)) .= $poly->((0))->clip($range->((0))->minmax);
    $poly->((1)) .= $poly->((1))->clip($range->((1))->minmax);
    $poly->((2)) .= $poly->((2))->clip($range->((2))->minmax);
    
    $ok = ($pol2 == $poly)->prodover;
    $rgb *= $ok->(*3);
    $prgb *= $ok->(*3);
    
    print "Calling points3d...\n" if($Flux::debug);
    points3d($range,zeroes($range),{PointSize=>0});
    hold3d;
    
    my $w3d = PDL::Graphics::TriD::get_current_window();

    my $nc = ($poly - $boxmin)/($boxmax-$boxmin);

    if($opt->{points} || !defined($opt->{points})) {
      my $p = new PDL::Graphics::TriD::Points($nc,$prgb,{PointSize=>($opt->{psize}||4)});
      $w3d->add_object($p);
    }

    ##############################
    # Generate line strips for each fluxon...
    my $rgbdex = 0;
    my @id = $w->fluxon_ids;
    print "Defining line strips...\n" if($Flux::debug);
    for my $i(0..$#id) {
	my $fp = $w->fluxon($id[$i])->polyline;
	my $normal_coords = ($fp-$boxmin)/($boxmax-$boxmin);

	my $l = new PDL::Graphics::TriD::LineStrip($normal_coords,$rgb[$i],{LineWidth=>($opt->{linewidth}||1)});
	$w3d->add_object($l);

#	line3d($fp,$rgb[$i]);

	if($opt->{label}) {
	    ##############################
	    # Label each point in the fluxon...
	    # This is really cheesy since label seems to take normalized coordinates 
	    # and I can't seem to figure the interface for normalizing plot coordinates.
	    # The normalization is done by start_scale, add_scale, and finish_scale in
	    # Graphics/TriD/Graph.pm within PDL, but they are somewhat opaque to me.
	    use PDL::Graphics::TriD::Labels;
	    my @labels = map { $_->id } $w->fluxon($id[$i])->vertices ;
		    
	    for my $j(0..$#labels) {
		my $l = new PDL::Graphics::TriD::Labels($normal_coords->(:,($j)),
							{Strings=>["$labels[$j]"],
							 Font=>$PDL::Graphics::TriD::GL::fontbase});
		$w3d->add_object($l);
	    }
	}
    }

    print "Neighbors?\n" if($Flux::debug);
    $nscale = $opt->{'nscale'} || 0.25;

    if($opt->{'neighbors'}){
	my @neighbors;

	for my $v (map { $_->vertices } $w->fluxons) {
	    next unless($v->next);
	    $xcen = 0.5 * ($v->x + $v->next->x);

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
	my $normal_coords = ($fp-$boxmin)/($boxmax-$boxmin);
	my $p = new PDL::Graphics::TriD::Points($normal_coords,{PointSize=>2});
	$w3d->add_object($p);
#	points3d(cat(@neighbors),{PointSize=>2});

    }

    if($opt->{'hull'}) {

      print "hullrgb...\n";
	$hullrgb = defined($opt->{'hullrgb'}) ? $opt->{'hullrgb'} : pdl(0.3,0.3,0);

      print "hullopen...\n";
	$hullopen = defined($opt->{'hullopen'}) ? $opt->{'hullopen'} : 10;

	for my $v( $w->vertices ) { ### map { $_->vertices } $w->fluxons) {
	    next unless($v->next);

	    $xcen = 0.5 * ($v->x + $v->next->x);
	    
	    my $pm = $v->projmatrix;
	    
	    my $hull = $v->hull;

	    my @hpoints = ();

	    for my $i(0..$hull->dim(1)-1) {

		my $rows = $hull->range([0,$i],[7,2],'p');
		if( ! ($rows->((4),:)->any) ) {
		    my $x = zeroes(3);
		    $x->(0:1) .= $rows->(2:3,(0));
		    $x *= $nscale;
		    $xx = $x x $pm;
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

			my $normal_coords = ($fp-$boxmin)/($boxmax-$boxmin);
			my $l = new PDL::Graphics::TriD::LineStrip($normal_coords,$hullrgb->dummy(1,$fp->dim(1))->copy);

			$w3d->add_object($l);

			@hpoints = ();
		    }


		    my $x = zeroes(3);
		    $x->(0) .= $rows->((2),(1)) + $hullopen * cos($rows->((6),(0)));
		    $x->(1) .= $rows->((3),(1)) + $hullopen * sin($rows->((6),(0)));
		    $x *= $nscale;
		    $xx = $x x $pm;
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
		    $xx = $x x $pm;
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

		hold3d;
		my $fp = cat(@hpoints)->(:,(0),:);
		my $normal_coords = ($fp-$boxmin)/($boxmax-$boxmin);


		my $l = new PDL::Graphics::TriD::LineStrip($normal_coords,$hullrgb->dummy(1,$fp->dim(1))->copy);
		$w3d->add_object($l);
	    }
	}    
    }



    print "ok.  twiddling...\n" if($twiddle);

    release3d;
    keeptwiddling3d  if($twiddle);
    print "twiddling...\n";
    twiddle3d();
    print "ok\n";
    nokeeptwiddling3d;

#    print "returning...\n";

}

=pod


=head2 polyline

=for ref

Dumps a giant polyline for all fluxons in the simulation.

=cut

sub polyline {
    my $w = shift;
    my @l,$l1;
    @l = map { $_->polyline } $w->fluxons;
    $l1 = @l[0];
    return $l1->glue(1,@l[1..$#l]);
}



=pod

=head2 bfield

Dumps the magnetic field value for all fluxons in the simulation

=cut

sub bfield {
    my $w = shift;
    my @l,$l1;
    @l = map { $_->bfield } $w->fluxons;
    $l1 = @l[0];
    return $l1->glue(1,@l[1..$#l]);
}


=pod

=head2 dump_vecs

=for ref

Dumps all vectors for all fluxons (see Flux::Fluxon::dump_vecs)

=cut

sub dump_vecs {
    my $w = shift;
    my @l,$l1;
    @l = map { $_->dump_vecs } $w->fluxons;
    $l1 = @l[0];
    return $l1->glue(1,@l[1..$#l]);
}

=pod

=head2 _hull_points

=for usage

    $hull = PDL::World::_hull_points($points);

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
    PDL::World::_plot_hull($w, PDL::World::_hull_points($points) [, $points]);

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
	for $i(0..$points->dim(1)-1){
	    $win->text($i,$points->((0),($i)),$points->((1),($i)),{color=>2});
	}
    }

    # Plot all points in neighbors.
    $win->points($hull->((0)),$hull->((1)),{color=>3});
    $win->hold;
    for $i(0..$hull->dim(1)-1){
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


sub DESTROY {
  # DESTROY gets called twice -- once for the tied hash and once for the underlying object
  # (which is a scalar ref, not a hash ref). Ignore the destruction for the tied hash.

  eval '$_[0]->{foo}';
  return unless($@);
  undef $@;

  _dec_refct_destroy( $_[0] );
}

######################################################################
# TIED INTERFACE
# Mostly relies on the general utility functions in Flux....

sub TIEHASH {
    my $class = shift;
    my $ptr = shift;
    my $me = \$ptr;
    return bless($me,$class);
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

