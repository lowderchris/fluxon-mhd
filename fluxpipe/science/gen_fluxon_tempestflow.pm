=head1 NAME

gen_fluxon_wsaflow - Comprehensive Module for Solar Wind Solution and Fluxon Flow Analysis

=cut

package gen_fluxon_tempestflow;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(gen_fluxon_tempestflow);

use PDL::Graphics::Gnuplot;
use PDL;
use PDL::NiceSlice;
use PDL::Options;
use PDL::ImageND;
use PDL ('which_min'); # This imports which_min directly into your current package
use POSIX;
use Math::RungeKutta;
use PDL::GSL::INTEG;
use PDL qw(squeeze);
use PDL::GSL::INTERP;

$PDL::verbose = 0;

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 USAGE


=head1 DEPENDENCIES

This module depends on the following Perl modules:
=over 4

=back

=head1 DIAGNOSTICS

Set the C<$verbose> flag to 1 for diagnostic output.

=head1 AUTHOR

Gilly <gilly@swri.org> (and others!)

=head1 SEE ALSO

Other relevant modules and documentation.

=cut



my $verbose = 0;

# Interpolate an array given X and Y values along with a test point
sub interpolate_1d {
    my ($x, $y, $x0) = @_;
    my $y0;
    my $len = $x->nelem;
    my $ind = abs($x - $x0)->minimum_ind;
    if ($ind == 0) {
        $y0 = $y->(($ind));
    } elsif ($ind == $len-1) {
        $y0 = $y->(($ind));
    } else {
        my $x1 = $x->(($ind));
        my $x2 = $x->(($ind+1));
        my $y1 = $y->(($ind));
        my $y2 = $y->(($ind+1));
        $y0 = $y1 + ($y2 - $y1) * ($x0 - $x1) / ($x2 - $x1);
    }
    return $y0;
}

sub gen_fluxon_tempestflow {
    my $fluxon = shift;
    my $fid = shift;
    my $output_file_name = shift;

    # Check if the start and end points of the fluxon are open
    my $is_start_open = ($fluxon->{fc_start}->{label} == -1);
    my $is_end_open = ($fluxon->{fc_end}->{label} == -2);

    # If the fluxon is labeled as a plasmoid, or if both ends are open, return undefined
    return undef if($fluxon->{plasmoid} || ($is_start_open + $is_end_open != 1));

    use File::Basename;

    # Get the file's base name and directory
    my ($name, $dir, $ext) = fileparse($output_file_name, qr/\.[^.]*/);

    # Replace .dat with _reinterpolated.dat
    my $file_path = "${dir}${name}_reinterpolated.dat";

    # Print the new filename for verification
    # print "$file_path\n";



    # print "\n\n\n\n\n";
    # print $output_file_name;
    # print "\n";
    # print "/Users/cgilbert/vscode/fluxons/fluxon-data/batches/tempest/data/cr2200/wind/cr2200_f1000_radial_bmag_tempest_reinterpolated.dat";
    # print "\n";
    # print $file_path;
    # print "\n\n\n\n\n";

    my $model_id = $fid;  # Change to desired model number
    my ($radii_ref, $velocities_ref) = get_tempest_velocity($file_path, $model_id);

    # # Output the results

    my $r_vr_scaled = pdl($radii_ref, $velocities_ref);

    (my $r_fr_scaled, my $theta, my $phis) =
            gen_fluxon_tempestflow_physics($fluxon, $fid);

    # Return the final fluxon array, fluxon radius, and magnetic field components
    return ($r_vr_scaled, $r_fr_scaled, $theta, $phis);
}




# Cache for file data to avoid repeated file reads
my %data_cache;

sub load_data {
    my ($filename) = @_;
    open my $fh, '<', $filename or die "Could not open file '$filename': $!";

    # Skip the header line
    my $header = <$fh>;

    # Read data
    while (my $line = <$fh>) {
        chomp $line;
        my ($fnum, $radius, $velocity) = split(/\s+/, $line);

        # Group data by model number (fnum)
        push @{$data_cache{$fnum}}, { radius => $radius, velocity => $velocity };
    }

    close $fh;
}

sub get_tempest_velocity {
    my ($filename, $model_number) = @_;

    # Load data from file if it hasn't been loaded yet
    load_data($filename) unless %data_cache;

    # Check if the model number exists and return arrays directly
    if (exists $data_cache{$model_number}) {
        my @radii = map { $_->{radius} } @{$data_cache{$model_number}};
        my @velocities = map { $_->{velocity} } @{$data_cache{$model_number}};
        return (\@radii, \@velocities);
    } else {
        warn "Model number $model_number not found in the data.";
        return;
    }
}


























# Function to read data from ZEPHYR file
sub read_data {
    my ($filename) = @_;

    # Open the file or die with an error message
    open my $fh, '<', $filename or die "Cannot open file $filename: $!";

    # Skip the header lines
    for (1..10) {
        <$fh>;
    }

    # Initialize arrays to store the data
    my (@r_Rsun, @rho, @u, @Valf, @T);

    # Read the file line by line
    while (my $line = <$fh>) {
        # chomp $line;
        my ($dummy, $r_Rsun_val, $rho_val, $u_val, $Valf_val, $T_val) = split /\s+/, $line;

        # # Store the data in arrays
        push @r_Rsun, $r_Rsun_val;
        push @rho, $rho_val;
        push @u, $u_val;
        push @Valf, $Valf_val;
        push @T, $T_val;
    }

    # Close the file handle
    close $fh;

    # Convert the arrays to PDL objects for numerical operations
    my $r_Rsun_pdl = pdl(\@r_Rsun);
    my $rho_pdl = pdl(\@rho);
    my $u_pdl = pdl(\@u) / 1e5;
    my $Valf_pdl = pdl(\@Valf) / 1e5;
    my $T_pdl = pdl(\@T);

    # Return the PDL objects
    return ($r_Rsun_pdl, $rho_pdl, $u_pdl, $Valf_pdl, $T_pdl);
}

sub interpolate_1d_array {
    my ($x, $y, $x_new) = @_;
    my $interp = PDL::GSL::INTERP->new("linear", $x, $y);
    my $y_new = $interp->interp($x_new);
    return $y_new;
}

sub gen_fluxon_tempestflow_physics {
    my $me = shift;
    my $fid = shift;
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

    # Pass options into local variables
    my $g0 = $opt{g0};
    my $r0 = $opt{r0};
    my $cs = $opt{cs} * 1e3; ## convert to m/s from km/s

    # Check for start and ending points being open
    our $en_open = ($me->{fc_end}->{label}==-2);

    # Calculate array of sphiserical coordinate positions and areas along the fluxon
    # Work along the correct direction depending on which end is open
    if($en_open) {
        my $x = squeeze($me->dump_vecs->(0));
        my $y = squeeze($me->dump_vecs->(1));
        my $z = squeeze($me->dump_vecs->(2));
        our $r1 = ($x**2 + $y**2 + $z**2)->sqrt * $opt{r0};
        our $r = 0.5 * $r1->range([[0],[1]],[$r1->dim(0)],'e')->sumover;
        our $thetas = acos($z/$r1*$opt{r0});
        our $phis = atan2($y, $x);
        our $A = pdl(map {$_->{A}} ($me->vertices));
        our $bfield = pdl(map {$_->{b_vec}} ($me->vertices));
        our $bmag = pdl(map {$_->{b_mag}} ($me->vertices));
        our $T = pdl(map {$_->{T}} ($me->vertices));
        our $rho = pdl(map {$_->{rho}} ($me->vertices));

    } else {
        my $x = squeeze($me->dump_vecs->(0,-1:0:-1));
        my $y = squeeze($me->dump_vecs->(1,-1:0:-1));
        my $z = squeeze($me->dump_vecs->(2,-1:0:-1));
        our $r1 = ($x**2 + $y**2 + $z**2)->sqrt * $opt{r0};
        our $r = 0.5 * $r1->range([[0],[1]],[$r1->dim(0)],'e')->sumover;
        our $thetas = acos($z/$r1*$opt{r0});
        our $phis = atan2($y, $x);
        our $A = pdl(map { $_->{A}} ($me->vertices))->(-1:0:-1);
        our $bfield = pdl(map {$_->{b_vec}} ($me->vertices))->(-1:0:-1);
        our $bmag = pdl(map {$_->{b_mag}} ($me->vertices))->(-1:0:-1);
        our $T = pdl(map {$_->{T}} ($me->vertices))->(-1:0:-1);
        our $rho = pdl(map {$_->{rho}} ($me->vertices))->(-1:0:-1);

    }

    our $A;
    our $bfield;
    our $bmag;
    our $r1;
    our $r;
    our $thetas;
    our $phis;
    our $T;
    our $rho;

    my $rn = $r1/$opt{r0};

    # Fix the bottom r vertex
    my $diff0 = abs(1 - $rn->((0)));
    $rn(0) .= 1 + $diff0;

    my $zn = $rn - 1;
    my $len = $rn->nelem - 1;


    # Get rid of end anomalies
    # $A->((0)) .= $A->((1));
    $A->((-1)) .= $A->((-2));

    # Fix the bottom magnetic field vertex
    $bmag->((1)) .= $bmag->((2)) * $A->((2)) / $A->((1));
    $bmag->((0)) .= $bmag->((1)) * $A->((1)) / $A->((0));

    ## Declare Variables
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



    # FIND THE SUN SURFACE
    my $first_ind = 1;
    my $r_sun = $rn->(($first_ind));
    my $B_sun = $bmag->(($first_ind));
    my $A_sun = $A->(($first_ind));

    my $phisi0 = $phis->(($first_ind));
    my $sin_theta0 = sin($thetas)->(($first_ind));
    my $r00 = $rn->(($first_ind));

    # FIND THE SOURCE SURFACE
    my $r_ss = 2.5; # Solar radii
    my $B_ss = interpolate_1d($rn, $bmag, $r_ss);
    my $A_ss = interpolate_1d($rn, $A, $r_ss);
    my $f_ss = abs($A_ss / $A_sun) * ($r_sun * $r_sun) / ($r_ss * $r_ss);

    # my $f_all= abs($B_sun / $bmag) * ($r_sun * $r_sun) /   ($rn * $rn);
    my $f_all =    abs($A / $A_sun) * ($r_sun * $r_sun) /   ($rn * $rn);

    # FIND THE TOP OF THE DOMAIN
    my $last_ind    = $phis->nelem - 1;
    my $phisi1      = $phis->(($last_ind));
    my $r11         = $rn->(($last_ind));
    my $sin_theta1  = sin($thetas)->(($last_ind));

    if ($phisi0 < 0) {$phisi0 += 2 * 3.1415926;}
    if ($phisi1 < 0) {$phisi1 += 2 * 3.1415926;}


    use strict;
    use warnings;
    use PDL;
    use PDL::Graphics::Gnuplot;

    # Sample data
    my $ones = ones($r->dims);


    # if (0) {
    #     # # Plot the data with labels
    #     my $plot = PDL::Graphics::Gnuplot->new(persist => 1);
    #     $plot->plot({title=>"Fluxon #$fid", logscale=>'xy', xlabel=>'z = r/R - 1', font=>'48'}, {legend=>"|B|"},$zn, $bmag);
    #     $plot->replot({legend=>"|B|", with =>'points'},$zn, $bmag);

    #     # $plot->replot({legend=>"bmag"},$zn, $bmag);
    #     # $plot->replot({legend=>"f_ss"},$rn, $ones * $f_ss);
    #     # $plot->replot({legend=>"f_all"} ,$rn, $f_all);
    #     $plot->replot({legend=>"f_r"},$zn, $f_all);
    #     $plot->replot({legend=>"f_r", with=>"points"},$zn, $f_all);

    #     $plot->replot({legend=>"A"},$zn, $A);
    #     $plot->replot({legend=>"A", with=>"points"},$zn, $A);

    #     $plot->replot({legend=>"d_f"},$zn, $densfac);
    #     $plot->replot({legend=>"d_f", with=>"points"},$zn, $densfac);

    #     $plot->replot({legend=>"u_w"},$zn, $u_scaled);
    #     $plot->replot({legend=>"u_w", with=>"points"},$zn, $u_scaled);

    #     # $plot->replot({legend=>"flow", with=>"points"},$rn->((-1)), $flow_field->((-1)));
    #     # $plot->replot({legend=>"phis", with=>"points"},$zn, $phis);

    #     # $plot->replot(legend=>"phis0", with=>"points", pointsize=>3, $r00, $phisi0);
    #     # $plot->replot(legend=>"phis1", with=>"points", pointsize=>3, $r11, $phisi1);

    #     # $plot->replot({legend=>"stheta"},$zn, sin($thetas));
    #     # $plot->replot({legend=>"stheta", with=>"points"}, $zn, sin($thetas));

    #     # $plot->replot(legend=>"sth0", with=>"points", pointsize=>3 ,$r00, $sin_theta0);
    #     # $plot->replot(legend=>"sth0", with=>"points", pointsize=>3 ,$r11, $sin_theta1);
    #     # $plot->replot({legend=>"distance"},$rn, $ones * $distance_degrees);

    #     # my $plot = PDL::Graphics::Gnuplot->new(persist => 1);
    #     # my $top_speed = $flow_field->((-1));
    #     # $plot->plot({title=>"ID: $fluxon_id, Top Speed: $top_speed", logscale => 'xy'}, {legend=>"rn", with=>"points"},$rn, $rn);

    #     # Block execution
    #     print "Press ENTER to continue...\n";
    #     <STDIN>;
    # }

    ## Calculate Return Values ##

    # our $r_v_scaled = pdl($rn, $u_scaled);

    my $r_fr_scaled = pdl($rn, $f_all);


    # Return the constructed arrays
    # return ($r_v_scaled, $r_fr_scaled, $thetas, $phis);
    return ($r_fr_scaled, $thetas, $phis);
}

1;
__END__
