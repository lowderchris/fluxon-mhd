=head2 make_world_sphere - given a list of lines, make a world.

=for ref

You feed in a list of pdls, one pdl per fluxon.  You get back a world
with default characteristics and two flux concentrations per fluxon (north
at the beginning, south at the end). All concentrations are on the spherical photosphere that you provide.


=cut

package make_world_sphere;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(make_world_sphere);
use PDL;
use PDL::NiceSlice;


sub make_world_sphere {
  my @lines = @_;
  my $opt = pop // {};
  push(@lines,$opt) and undef $opt unless UNIVERSAL::isa($opt,'HASH');
  ##push it back onto @lines if it isn't a hash

  my $phot=$opt->{photosphere};
  my $rad=$opt->{rad} // 1;
  print "radius is $rad \n";

  @lines = @_;
  push(@lines,$opt) and undef $opt unless UNIVERSAL::isa($opt,'HASH');
  ##push it back onto @lines if it isn't a hash

  my @out = "";

  push(@out, "GLOBAL FORCES f_pressure_equi2b f_curvature f_vertex4");
  push(@out, "GLOBAL BOUNDARY SPHERE 0 0 0 ".$rad);
  push(@out, "GLOBAL OPEN 0 0 0 ".($opt->{rmax}+1)." 1") if($opt->{rmax});
  push(@out,"");

  my $j=100;
  my $ln=0;

  foreach my $l(@lines) {
      my $line=$l->copy;
      my ($fc0, $fc1, $fl);

      my $rstart = sqrt(sum($line->(:,(0)) * $line->(:,(0))));
      my $rend   = sqrt(sum($line->(:,(-1)) * $line->(:,(-1))));

      ##simply scale the first/last to the right distance. assume that we have the right # of points.
      $line->(:,(0)) *= $rad / sqrt(sum($line->(:,(0))*$line->(:,(0))));
      $line->(:,(-1)) *= $rad / sqrt(sum($line->(:,(-1))*$line->(:,(-1))));

##      if($rstart > $rad) {
##	  $line->(:,(0)) *= $rad / sqrt(sum($line->(:,(0))*$line->(:,(0))));
##      } else {
##	  ##assume that only that last one is under the line and do the simple thing, scale the last poiont.
##	  if (sqrt(sum($line->(:,(1))*$line->(:,(1)))) < $rad) {
##	      print "uh-oh, more than one point is under the sphere";
##	  }

	  ##my $a=sum($line->(:,(1))*$line->(:,(1)));
	  ##my $b=sum( ($line->(:,(0))-$line->(:,(1)))*($line->(:,(0))-$line->(:,(1))) );
	  ##my $c=sum( $line->(:,(1))*($line->(:,(0))-$line->(:,(1))) );

##      if($rstart < $rad) {
##	  my $rn = sqrt(sum($line->(:,(1))*$line->(:,(1))));##next rad
##	  if($rn > $rad) {
##	      my $alpha = ($rad - $rstart) / ($rn - $rstart);
##	      my $newv = $line->(:,(0)) * ($rad-$alpha) + $line->(:,(1)) * $alpha;
##	      $line->(:,(0)) .= $newv;
##	  }
##	  # above uses planar photosphere, not quite right.   Fudge the last bit.
##	  $line->(:,(0)) /= sqrt(sum($line->(:,(0))*$line->(:,(0))));
##      }

##    if($rend < 1) {
##	my $rp = sqrt(sum($line->(:,(-2))*$line->(:,(-2))));
##	if($rp > 1) {
##	    my $alpha = (1 - $rend) / ($rp - $rend);
##	    my $newv = $line->(:,(-1)) * (1-$alpha) + $line->(:,(-2)) * $alpha;
##	    $line->(:,(-1)) .= $newv;
##	}
##	$line->(:,(-1)) /= sqrt(sum($line->(:,(-1))*$line->(:,(-1))));
##    }

    my $open_start = ($opt->{rmax} && $rstart >= $opt->{rmax});
    my $open_end   = ($opt->{rmax} && $rend   >= $opt->{rmax});
      ##p $open_start, $open_end, "\n";

    unless($open_start){
      push(@out,sprintf("NEW %d %9.3g %9.3g %9.3g 1",
			($fc0=$j++),
			$line->at(0,0), $line->at(1,0), $line->at(2,0)
			)
	   );
    }

    unless($open_end){
      push(@out,sprintf("NEW %d %9.3g %9.3g %9.3g -1",
			($fc1=$j++),
			$line->at(0,-1), $line->at(1,-1), $line->at(2,-1)
			)
	   );
    }

    push(@out,sprintf("LINE %d %d %d 1",
		      ($fl = $j++),
		      ($open_start ? -1 : $fc0),
		      ($open_end ? -2 : $fc1)
		      )
	 );
    $ln++;
    for my $k(1..$line->dim(1)-2) {
     ## print "line $ln; k=$k\n";
      push(@out,sprintf("VERTEX %d %d %d %9.3g %9.3g %9.3g",
			$fl,
			$j++,
			$k,
			$line->(:,($k))->list
			)
	   );
    }
    push(@out,"");
  }
  join("\n",@out);
}
1;