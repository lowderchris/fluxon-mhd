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
package Flux::World;

require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader);
@EXPORT = qw( read_world write_world);

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
    return "Fluxon geometry object; ".(@lines+0)." fluxons, ".(@vertices+0)." vertices\n";
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
	    _set_force($me,$i,$s); # Set forces one at a time (in World.xs)
	}
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
  
1;



