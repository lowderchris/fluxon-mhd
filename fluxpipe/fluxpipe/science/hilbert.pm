=head2 hilbert - generate a Hilbert ordering of the points in a PDL

=for usage

$line = hilbert($pdl);

$line = hilbert($w,$h);

=for ref

The Hilbert curve fills the unit square with subsequent refinements,
and yields a pixel ordering of the plane that is better than rastering
for some purposes.  The output of hilbert is a 2xN PDL that contains
coordinates of all the points in $pdl, or in the square (0..$x-1,
0..$y-1), in the order that they are visited by a Hilbert curve of
order ceil(log2(max($w,$h))).

The nth Hilbert curve is the nth element in a convegent series of curves
used to produce a fractal with Hausdorff dimension 2.  The refined curves
are crinkly and hence have two benefits:  they are (more or less) local, in
the sense that points with close indices on the Hilbert curve tend to be
close in the plane; and they are pseudorandom in the sense that it is
difficult to predict the spatial offset between two nonadjacent points
on the curve.

Since Hilbert curves expand only to 2^N pixels, the algorithm chooses the
next larger 2^N x 2^N square, and then ignores the points that happen to
be outside the area you wanted.  That yields a nonlocal ordering of the
points in the rectangle, because some points (near the boundaries) aren't
close (in 2-D) to their neighbors in the ordering.

If you want the curve to be as local as possible, make sure that w and h
are both the same power of 2.

hilbert uses the Inline::Pdlpp package to do the heavy lifting in C.

=cut

package hilbert;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(hilbert);
use PDL;
use PDL::NiceSlice;

sub hilbert {

  my ($w, $h, $keep);
  if(UNIVERSAL::isa($_[0], 'PDL')) {
      $w=$_[0]->dim(0);
      $h=$_[1]->dim(1);
  } else {
      $w = shift;
      $h = shift;
  }

  my $dim = ( $w>$h ? $w : $h );

  my $n = (log(pdl($dim))/log(2))->ceil;
  my $siz = 2 ** $n;

  my $coords = zeroes(long, 2, $siz*$siz);

  PDL::hpp($coords,$n);

  my $of = pdl(long, ($siz-$w)/2, ($siz-$h)/2);

  if($keep) {
      return $coords - $of;
  } else {
      $coords -= $of;
      my $c2w = which( ($coords->((0)) < $w) & ($coords->((1)) < $h) &
		    ($coords->((0)) >= 0) & ($coords->((1)) >= 0));
      return $coords->(:,$c2w)->sever;
  }
}

no PDL::NiceSlice;
use Inline Pdlpp=><<'EOF';

pp_addhdr( << 'EOAHD');
#define H_UP    '^'
#define H_LEFT  '<'
#define H_RIGHT '>'
#define H_DOWN  'v'
static char *hilbert_recurse(char *where, char dir, int level) {
    if(level==0) {
	switch(dir) {
	    case H_LEFT:
	    *(where++) = H_RIGHT;
	    *(where++) = H_DOWN;
	    *(where++) = H_LEFT;
	    break;

	    case H_RIGHT:
	    *(where++) = H_LEFT;
	    *(where++) = H_UP;
	    *(where++) = H_RIGHT;
	    break;

	    case H_UP:
	    *(where++) = H_DOWN;
	    *(where++) = H_RIGHT;
	    *(where++) = H_UP;
	    break;

	    case H_DOWN:
	    *(where++) = H_UP;
	    *(where++) = H_LEFT;
	    *(where++) = H_DOWN;
	    break;
	}
    } else {
	int ln = level - 1;
	switch(dir) {
	    case H_LEFT:
	    where = hilbert_recurse(where, H_UP, ln);
	    *(where++) = H_RIGHT;
	    where = hilbert_recurse(where, H_LEFT, ln);
	    *(where++) = H_DOWN;
	    where = hilbert_recurse(where, H_LEFT, ln);
	    *(where++) = H_LEFT;
	    where = hilbert_recurse(where, H_DOWN, ln);
	    break;

	    case H_RIGHT:
	    where = hilbert_recurse(where, H_DOWN, ln);
	    *(where++) = H_LEFT;
	    where = hilbert_recurse(where, H_RIGHT, ln);
	    *(where++) = H_UP;
	    where = hilbert_recurse(where, H_RIGHT, ln);
	    *(where++) = H_RIGHT;
	    where = hilbert_recurse(where, H_UP, ln);
	    break;

	    case H_UP:
	    where = hilbert_recurse(where, H_LEFT, ln);
	    *(where++) = H_DOWN;
	    where = hilbert_recurse(where, H_UP, ln);
	    *(where++) = H_RIGHT;
	    where = hilbert_recurse(where, H_UP, ln);
	    *(where++) = H_UP;
	    where = hilbert_recurse(where, H_RIGHT,ln);
	    break;

	    case H_DOWN:
	    where = hilbert_recurse(where, H_RIGHT, ln);
	    *(where++) = H_UP;
	    where = hilbert_recurse(where, H_DOWN, ln);
	    *(where++) = H_LEFT;
	    where = hilbert_recurse(where, H_DOWN, ln);
	    *(where++) = H_DOWN;
	    where = hilbert_recurse(where, H_LEFT, ln);
	    break;
	}
    }
    return where;
}
EOAHD

pp_def('hpp',
       Pars=>'out(n,m);',
       GenericTypes=>[L],
       OtherPars=>'int level',
       Code=> <<'EOHPP'
       long dex;
       long side;
       long x;
       long y;
       char *string = (char *)malloc($SIZE(m)+2);
       char *s = string;

       s = hilbert_recurse(string, H_RIGHT, $COMP(level)-1);
       *(s++) = 0;

       side = 1 << $COMP(level);
       s=string;

       loop(m) %{
	   if(!m){
	       x = 0;
	       y = 0;
	   } else {
	       switch( *(s++) ){
		 case H_UP:    y++; break;
		 case H_DOWN:  y--; break;
		 case H_LEFT:  x++; break;
		 case H_RIGHT: x--; break;
	       }
	       $out( n=>0 ) = x;
	       $out( n=>1 ) = y;
	   }
       %}

       free(string);

EOHPP
       );



EOF
