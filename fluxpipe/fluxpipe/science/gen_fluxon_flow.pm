use PDL::NiceSlice;

=head1 NAME

gen_fluxon_flow - Calculate Flow and Density in a Fluxon

=cut

package gen_fluxon_flow;
use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw(gen_fluxon_flow);

use PDL::Options;  ## This is where the 'parse' function comes from
use PDL::ImageND;  ## This is where the ConvolveND comes from
use Math::RungeKutta; ## This requires extra package: cpan install Math::RungeKutta
use PDL::GSL::INTEG; ## This requires extra package: cpan install PDL::GSL::INTEG, and ALSO you have to install GSL to the machine, first.
use PDL qw(squeeze);

=head1 SYNOPSIS

Given a fluxon and some options, this function calculates the flow and density in the fluxon if it is open,
or leaves it blank if it is closed.

=head1 DESCRIPTION

This function returns a 2xN PDL containing (r, vr). By default, it also updates these fields in the fluxon itself.
Additionally, it returns two PDLs containing the starting and ending theta and phi values for the fluxon.

If the fluxon is doubly open, doubly closed, or a plasmoid, the function returns undef.

=head1 USAGE

    use gen_fluxon_flow;
    my $result = gen_fluxon_flow($fluxon, \%options);

=head1 OPTIONS

=over 4

=item * steps: Number of steps for calculation. Default is 500.

=item * g0: Gravitational acceleration in m/s^2. Default is 280.

=item * r0: Radius in meters. Default is 696e6.

=item * v0: Velocity in km/sec. Default is a random number between 10 and 20.

=item * cs: Sound speed in km/sec at 2e6K and fully ionized protons. Default is 180.

=back

=head1 DEPENDENCIES

This module depends on the following Perl modules:

=over 4

=item * PDL::NiceSlice
=item * PDL::Options
=item * PDL::ImageND
=item * Math::RungeKutta
=item * PDL::GSL::INTEG

=back

To install the extra packages, you can use:

    cpan install Math::RungeKutta
    cpan install PDL::GSL::INTEG

You also need to install GSL on your machine.

=head1 AUTHOR

Gilly <gilly@swri.org> (and others!)


=cut


sub gen_fluxon_flow {
    my $me = shift;
    my $u_opt = shift // {};

    # Define and append optional option input
    my %opt = parse( {
        steps => 500,
        g0 => 280,      # m/s^2
        r0 => 696e6,    # m
        #v0 => 10,       # km/sec
        v0 => rand(10)+10,       # km/sec
        cs => 180,      # km/sec at 2e6K and fully ionized protons
                        },
                        $u_opt
        );

    # print $opt{v0};
    # Check for start and ending points being open
    our $st_open = ($me->{fc_start}->{label}==-1);
    our $en_open = ($me->{fc_end}->{label}==-2);

    # Return undefined value if the given fluxon structure is labeled as a plasmoid, or if
    #   both ends are open
    return undef if( $me->{plasmoid} || ($st_open + $en_open != 1));

    # Pass options into local variables
    my $g0 = $opt{g0};
    my $r0 = $opt{r0};
    my $cs = $opt{cs} * 1e3; ## convert to m/s from km/s

    # Calculate array of spherical coordinate positions and areas along the fluxon
    # Work along the correct direction depending on which end is open
    if($en_open) {
        my $x = squeeze($me->dump_vecs->(0));
        my $y = squeeze($me->dump_vecs->(1));
        my $z = squeeze($me->dump_vecs->(2));
        our $r1 = ($x**2 + $y**2 + $z**2)->sqrt * $opt{r0};
        our $r = 0.5 * $r1->range([[0],[1]],[$r1->dim(0)],'e')->sumover;
        our $th = acos($z/$r1*$opt{r0});
        our $ph = atan2($y, $x);
        our $A = pdl(map {$_->{A}} ($me->vertices)) * $opt{r0} * $opt{r0};
    } else {
        my $x = squeeze($me->dump_vecs->(0,-1:0:-1));
        my $y = squeeze($me->dump_vecs->(1,-1:0:-1));
        my $z = squeeze($me->dump_vecs->(2,-1:0:-1));
        our $r1 = ($x**2 + $y**2 + $z**2)->sqrt * $opt{r0};
        our $r = 0.5 * $r1->range([[0],[1]],[$r1->dim(0)],'e')->sumover;
        our $th = acos($z/$r1*$opt{r0});
        our $ph = atan2($y, $x);
        our $A = pdl(map { $_->{A}} ($me->vertices))->(-1:0:-1) * $opt{r0} * $opt{r0};
        $A(0:-2) .= $A(1:-1);
    }

    ## Declare Variables
    our $A;
    our $r;
    our $th;
    our $ph;
    # our $vu;
    my $vnow;
    my $dr;
    my @rv;
    my $rnow;
    my @num;
    my @denom;
    my $step;
    my $vint;
    my $abserr;
    my $ierr;

    # Get rid of end anomalies
    $A->((0)) .= $A->((1));
    $A->((-1)) .= $A->((-2));


    # Define gravitational acceleration g as a function of radius
    our $g = $g0 * ($opt{r0} / $r)**2;

    # Define an expansion factor
    # Note the direction has been set in the reading
    our $fr = $A / $A->((1)) * ($r->((1))) * ($r->((1))) / $r / $r;

    # Smooth slightly to remove numerical issues
    $fr->inplace->convolveND(pdl(1,2,3,2,1)/9,{bound=>'e'});
    our $dfr_dr = $fr->copy;

    # Calculate numerical derivatives
    $dfr_dr->(1:-2) .= ($fr->(2:-1) - $fr->(0:-3))/($r->(2:-1) - $r->(0:-3));
    $dfr_dr->(0) .= ($fr->(1) - $fr->(0)) / ($r->(1) - $fr->(0));
    $dfr_dr->(-1) .= ($fr->(-1) - $fr->(-2)) / ($r->(-1) - $r->(-2));

    # Chunk of code to calculate numerical derivative
    my $dydx = sub {
        my $rnow = shift;
        #my $vnow = shift;

        my $front = $vnow/$rnow;
        my $denominator = $cs * $cs  -  $vnow * $vnow;

        my ($gnow,$gerr) = interpolate($rnow,$r,$g);
        my ($frnow,$frerr) = interpolate($rnow,$r,$fr);
        my ($dfr_drnow,$dfrerr) = interpolate($rnow,$r,$dfr_dr);
        #my $dfr_drnow = 0;

        my $numerator = $gnow * $rnow   -   2 * $cs * $cs * (  1   +   $rnow / 2 / $frnow * $dfr_drnow  );
        #return ($front * $numerator / $denominator, $numerator, $denominator);
        return $front * $numerator / $denominator;
    };

    # Calculate average steps in r
    $dr = ($r->max - $r->min) / $opt{steps};

    # Define v0 and r0 values, and push into a list for storage
    @rv = ();
    $vnow = $opt{v0} * 1e3; # Convert to m/sec from km/sec
    $rnow = $r->minimum;
    push(@rv, pdl($rnow,$vnow));

    # Create some empty list for storage of numerators and denominators
    @num = ();
    @denom = ();

    my ($dvdr, $oldnum, $olddenom) = &$dydx($rnow,$vnow);

    for ($step=0; $step<$opt{steps};$step++){
        #($rnow, $vnow) = rk4([$vnow], \&$dydx, $rnow, $dr);
        ($vint, $abserr, $ierr) = gslinteg_qags($dydx, $rnow, $rnow+$dr, 0, 1e-6, 1000,{Warn => 'y'});
        $rnow += $dr;
        $vint = squeeze($vint);
        $vnow += $vint;
        push(@rv, pdl($rnow,$vnow));
    }

    ## Calculate Return Values ##

    # This array is (r, v) in units of (r_sun, km/s)
    my $r_v_scaled = pdl(@rv) * pdl(1.0 / $opt{r0}, 1e-3);

    # This array is (r, fr) in units of (r_sun, unitless)
    my $r_fr_scaled = pdl($r / $opt{r0}, $fr);

    # Return the constructed arrays
    return ($r_v_scaled, $r_fr_scaled, $th, $ph);

}
1;
__END__
