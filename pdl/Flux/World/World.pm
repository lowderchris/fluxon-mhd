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

=head1 FUNCTIONS

=cut

BEGIN {

use PDL::Graphics::TriD; 

package Flux::World;

require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader);
@EXPORT = qw( read_world write_world str2world );

bootstrap Flux::World;
}

package Flux::World;
use overload '""' => \&_stringify;



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

Produce a short summary string describing a fluxon model

This is the overloaded stringification routine for Flux::Worlds, so
using a World in string context gives you the short summary.  It is
implemented VERY cheesily and slowly (dumps the whole world to ASCII,
then checks how many LINE lines and VERTEX lines there are), so, er,
fix that.

=cut
    
sub _stringify {
    my $me = shift;
    my @s = $me->string;
    my @lines = grep(m/^\s*LINE/,@s);
    my @vertices = grep(m/^\s*VERTEX/,@s);
    return "Fluxon geometry object; ".(@lines+0)." fluxons, ".(@vertices+0)." vertices\n".
	"\tForces (".($me->_b_flag?"":"not")." B-normalized): ".join(", ",$me->forces)."\n";

}

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

    $baddies = $world->fix_curvature($thresh);

=for ref

Adds additional vertices as needed to ensure that no vertices have too 
much curvature.  The $thresh is the maximum bend angle, in radians, that
is allowed at each vertex.  You get back a count of the number of fluxels that
required subdivision.

You have to have primed the world with at least one timestep or
update_neighbors call before you do this -- otherwise horrible things might
happen.

=cut

# Implemented in World.xs



=pod

=head2 update_neighbors

=for usage

    $world->update_neighbors($globalflag);

=for ref

Update the neighbor trees between adjacent fluxels.  If $globalflag is
nonzero, then all previous neighbor information is discarded, causing
a global neighbor search for each fluxel.  The global algorithm is
O(n^2).  If $globalflag is zero, then a next-nearest-neighbors
algorithm is used, which is O(n).  The first time you call this, you 
had better set $globalflag=1.

=cut

# Implemented in World.xs


=pod

=head2 verbosity

=for usage

    $v = $world->verbosity();
    $world->verbosity($v);

Controls the level of, well, verbosity associated with operations on $world.
High values send insane amounts of stuff out the stderr port.

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

Return the fluxon(s) with the given IDs, as a list of Flux::Fluxon objects.

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

=head2 forces

=for usage
 
  print $world->forces;
  $world->forces('t_curvature','t_vertex')

=for ref

Return the names of the forces that are queued up to act on each fluxel
in the simulation, or set them.  The forces have predefined names
that are present at the top of physics.c.  Unknown forces are referred to
by hex value, but that way lies madness.

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
When set, the flag causes all the forces to be treated as if the 
local magnetic field magnitude were divided out.  This is appropriate
for the older forces (f_pressure_equi and f_curvature) but not
for the newer ones.

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

=head2 photosphere

=for usage

 $phot = pdl($world->photosphere);
 $world->photosphere([$phot->list]);

=for ref

With no additional arguments, returns a six-element list containing the 
(origin, normal-vector) coordinates of the photospheric plane, or a
zero-element list indicating that there is no photospheric plane.  If you
pass in a list ref, it is used to set the photospheric plane parameters.

=cut

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

Updates positions of the vertices after the forces have been calculated.
Takes one step toward relaxation.  Warning -- if $dt is too big you'll 
be in trouble, guv!

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

 $world->render_lines

=for ref

Produces a simple 3-D plot of all the field lines in a World, using PDL::TriD.
This bit of code has the potential to expand without limit -- it might be good to break
it out into a separate file....

=cut

use PDL::Graphics::TriD;
use PDL;
use PDL::NiceSlice;

sub render_lines {
    my $w = shift;
    my $twiddle = shift;
    my $range=shift;
    my $opt=shift;
    
    if(!defined $opt && ref $range eq 'HASH') {
	$opt = $range;
	undef $range;
    }


    release3d;

    my $fid;
    my @poly = map { $w->fluxon($_)->polyline } $w->fluxon_ids;

    my @rgb,@prgb;
    if($opt->{'RGB_CODE'}) {
	eval $opt->{'RGB_CODE'};
    } else {
	@rgb = map { ones($_) - (yvals($_)==0) } @poly;
    }
    if($opt->{'PRGB_CODE'}) {
	eval $opt->{'PRGB_CODE'};
    } else {
	@prgb = map { ones($_)*pdl(1,0,0) } @poly;
    }
    
    $poly = pdl(0,0,0)->glue(1,@poly);
    $rgb = pdl(0,0,0)->glue(1,@rgb);
    $prgb = pdl(0,0,0)->glue(1,@prgb);


    nokeeptwiddling3d;

    if(defined $range) {
	$pol2 = $poly->copy;
	$poly->((0)) .= $poly->((0))->clip($range->((0))->minmax);
	$poly->((1)) .= $poly->((1))->clip($range->((1))->minmax);
	$poly->((2)) .= $poly->((2))->clip($range->((2))->minmax);

	$ok = ($pol2 == $poly)->prodover;
	$rgb *= $ok->(*3);
	$prgb *= $ok->(*3);

      points3d($range,zeroes($range),{PointSize=>0});
      hold3d;
      line3d($poly,$rgb);
    } else {
      line3d($poly,$rgb);
      hold3d;
    }
    points3d($poly,$prgb,{PointSize=>4});
    release3d;

    keeptwiddling3d  if($twiddle);
    twiddle3d();
    nokeeptwiddling3d;
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
  
1;



