=head1 NAME

FluxSystem - Fluxon MHD simulation arena

=head1 SYNOPSIS

use PDL;
use FluxSystem;

$world = new FluxSystem;

$world = FluxSystem->ReadFromString($string);
$world = FluxSystem->ReadFromFile($filename);

$str = $world->WriteToString();
$ok  = $world->WriteToFile($filename);

$world->step($dt);
$lines = $world->vectors;

=head1 DESCRIPTION

FluxSystem objects are the perl interface to FLUX MHD models.  A $world
object is just a pointer off into hyperspace.  The perl destructor causes 
the fluxon to free itself.

=head1 FUNCTIONS

=cut

BEGIN {
package FluxSystem;

require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader);
@EXPORT = qw( read_world write_world);

bootstrap FluxSystem;
}

package FluxSystem;
use overload '""' => \&_stringify;

=head2 read_world

=for usage

$world = read_world($file);

=for ref

Read in a fluxon structure from a file

You supply the file name; the fluxon gets read into memory.
The file is a "standard" FLUX ASCII fluxon description file.

=head2 str2world

=for usage

$world = str2world($str)
$world = str2world(@lines)
$world = str2world([@lines])

=for ref

Convert a string or collection of lines to a fluxons structure

This allows you to convert a perl string to a fluxon structure.  It
is implemented cheesily:  the 'c' code uses streams to interpret the data,
so your strings are written to a temporary file that is then snarfed back
in.

=cut

sub str2world {
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
    
sub _stringify {
    my $me = shift;
    my @s = $me->string;
    my @lines = grep(m/^\s*LINE/,@s);
    my @vertices = grep(m/^\s*VERTEX/,@s);
    return "Fluxon geometry object; ".(@lines+0)." fluxons, ".(@vertices+0)." vertices\n";
}

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

1;



