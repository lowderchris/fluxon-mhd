=head1 NAME

Flux::Vertex - fluxon vertex / fluxel

=head1 SYNOPSIS

  use PDL;
  use Flux;
  $world = new Flux::World;
  
  $fluxon = $world->fluxon(19);
  $vertex = $fluxon->vertex(4);
  $vertices = $fluxon->vertices;

=head1 Methods

=cut

BEGIN {
package Flux::Vertex;

require Exporter;
require DynaLoader;
@ISA = qw( Exporter DynaLoader );
@EXPORT = qw( );

bootstrap Flux::Vertex;
}

package Flux::Vertex;

use overload '""' => \&stringify;

=pod

=head2 stringify

=for ref 

Produce a string summary of a vertex ("" overload)

=cut

sub stringify {
    &_stringify(shift);
}



=pod

=head2 id

=for ref

Return the integer id of the vertex

=cut

# Implemented in Vertex.xs

=pod

=head2 fluxon

=for ref

Return the fluxon associated with the vertex

=cut

# Implemented in Vertex.xs


=pod

=head2 prev

=for ref

Return the previous vertex on the same fluxon

=head2 next

=for ref

Return the next vertex on the same fluxon

=cut

# Implemented in Vertex.xs



=pod

=head2 neighbors

=for ref

Return the neighbors of this vertex, as a perl list of VERTEX objects.

=head2 nearby

=for ref

Return the inverse-neighbor list of this vertex, as a perl list of VERTEX objects.

=cut

sub neighbors {
  my $this = shift;
  return @{_adjacent($this,0)};     # in Vertex.xs
}

sub nearby {
  my $this = shift; 
  return @{_adjacent($this,1)};     # in Vertex.xs
}

=pod

=head2 hull

=for ref

Return the hull of a vertex, in 2-D.

Calculates the neighborhood and 2-D projected hull, and returns the 
coordinates as a PDL.  The output is a 5xN PDL.  The columns are:
[2d-x, 2d-y, hull-x, hull-y, open]
where 2d-x and 2d-y are the projected coordinates of the neighbor points,
and hull-x and hull-y are the coordinates of the hull vertices.
"open" is a flag indicating whether the associated vertex is "open" (ie
nonexistent)

=cut

# Implemented in Vertex.xs

=head2 proj_neighbors

=for usage

    $xyl = $vertex->proj_neighbors(0); # just neighbor vertices
    $xyl = $vertex->proj_neighbors(1); # all vertices in the World

=for ref

Return just the neighbors, projected into the perpendicular plane of the 
vertex's segment using the fancy projection.  The output is a 3xN PDL.
The columns are: [2d-x,2d-y,label] where 'label' is the label of the other vertex.

Takes a single parameter - a global flag to indicate that all vertices, not just the
precalculated neighbors, are to be included in the calculation.

=cut

# Implemented in Vertex.xs

=head2 plot_neighbors

=for usage

    $v->plot_neighbors($window, $global, $hull, [$env_opt]);

=for ref

Makes a 2-D plot of the neighbors of a particular vertex, with or without hull vertices.
The $window should be an open PGPLOT window.  This is useful mainly for debugging but is
a handy tool.  $global is a global flag to be handed to proj_neighbors; $hull indicates
whether the hull should be plotted on the same diagram; and $env_opt, if present, is an 
options hash to be passwd into the PGPLOT plotting routine.

=cut
use PDL::NiceSlice;

sub plot_neighbors {
    my($me, $window, $global, $hull, $env_opt) = @_;

    $window->release;

    my $xyl = $me->proj_neighbors($global);
    $Flux::xyl = $xyl;
    
    $env_opt = {} unless (defined $env_opt and ref $env_opt eq 'HASH');
    $env_opt->{title}="Projected neighbors of vertex ".$me->id 
	unless defined($env_opt->{title});

    $window->points($xyl->((0)),$xyl->((1)),$env_opt);
    $window->hold;

    for my $i(0..$xyl->dim(1)-1){
	$window->text($xyl->at(2,$i),$xyl->at(0,$i),$xyl->at(1,$i),{color=>3,charsize=>0.67,justification=>0.5*int(rand(3))});
    }

    if($hull) {
	my $hul = $me->hull();
	$Flux::hul = $hul;
	for $i(0..$hul->dim(1)-1) {
	    my $row = $hul->(:,($i));
	    my $r2 = $hul->range([0,$i+1],[7,0],'p')->sever;
	    if($row->((4))) {
		print "open before\n";
		# open vertex
		$x0 = cos($row->at(6))*1000;
		$y0 = sin($row->at(6))*1000;
		$window->line([$x0,$r2->at(2)],[$y0,$r2->at(3)],{color=>2});
	    } elsif($r2->((4))) {
		print "open after\n";
		$x0 = cos($r2->at(5))*1000;
		$y0 = sin($r2->at(5))*1000;
		$window->line([$row->at(2),$x0],[$row->at(3),$y0],{color=>2});
	    } else {
		$hv = $hul->range([2,$i],[2,2],'p');
		$window->line($hv->((0)),$hv->((1)),{color=>2});
	    }
	}
    }
}


=pod

=head2 projmatrix

=for ref

Returns the projection matrix to convert 3-vectors into the perpendicular plane to the segment.

=cut

# Implemented in Vertex.xs

=pod

=head2 x

=for ref

Returns the location vector of the vertex, as a 3-pdl.

=cut

# implemented in Vertex.xs

1;
