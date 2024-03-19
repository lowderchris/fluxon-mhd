=head1 NAME

gen_fluxon_wsaflow - Comprehensive Module for Solar Wind Solution and Fluxon Flow Analysis

=cut

package gen_fluxon_schonflow;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(gen_fluxon_schonflow do_image_plot);

use PDL;
use PDL::NiceSlice;
use PDL::Options;
use PDL::ImageND;
use Math::RungeKutta;
use PDL::GSL::INTEG;
use PDL qw(squeeze);
$PDL::verbose = 0;

=head1 SYNOPSIS

The module provides two primary functionalities:
1. gen_fluxon_wsaflow: Iteratively generates a transonic solar wind solution for a given fluxon.
2. gen_fluxon_wsaflow_phys: Calculates the flow and density in a fluxon, depending on its open or closed status.

=head1 DESCRIPTION

gen_fluxon_wsaflow - Takes a fluxon and optional user-defined settings as input, returning PDL arrays representing the final fluxon array, fluxon radius, and magnetic field components. Handles special cases like doubly open, doubly closed, or plasmoids.

gen_fluxon_wsaflow_phys - Provides detailed flow and density calculations in a fluxon, returning a 2xN PDL containing (r, vr), and updates these fields in the fluxon. Also returns theta and phi values for the fluxon, handling special cases similarly.

=head1 USAGE

For gen_fluxon_wsaflow:
    use gen_fluxon_wsaflow;
    my ($r_vr_scaled, $r_fr_scaled, $theta, $phi) = gen_fluxon_wsaflow($fluxon, \%options);

For gen_fluxon_wsaflow_phys:
    use gen_fluxon_wsaflow_phys;
    my $result = gen_fluxon_wsaflow_phys($fluxon, \%options);

=head1 OPTIONS

gen_fluxon_wsaflow options:
=over 4
=item * initial_temperature: Initial temperature in Kelvin. Default is 1e6.
=item * sound_speed: Speed of sound based on the initial temperature. Calculated internally.
=item * velocity_increment: Incremental step for velocity search. Default is 25 km/s.
=back

gen_fluxon_wsaflow_phys options:
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

=head1 DIAGNOSTICS

Set the C<$verbose> flag to 1 for diagnostic output.

=head1 AUTHOR

Gilly <gilly@swri.org> (and others!)

=head1 SEE ALSO

Other relevant modules and documentation.

=cut



my $verbose = 0;

sub interpolate_2d_lonlat {
    our ($image, $long_i, $latt_i) = @_;

    $image = pdl($image);
    # Define your image dimensions
    my ($img_width, $img_height) = $image->dims; # 900 by 360

    # Create grids for latitude and longitude
    my $long_vals = pdl(sequence($img_width)/($img_width-1) * 2 * 3.14159265);  # 0 to 2*pi
    my $latt_vals = pdl(sequence($img_height)/($img_height-1) * 2 - 1);  # -1 to 1

    my ($ind_long) = minimum_ind(abs($long_vals - $long_i));
    my ($ind_latt) = minimum_ind(abs($latt_vals - $latt_i));
    my $imval = $image->at($ind_long, $ind_latt);
    return $imval;
}

sub gen_fluxon_schonflow {
    my $fluxon = shift;
    my $distance_array_degrees = shift;
    my $fluxon_id = shift;

    # Check if the start and end    points of the fluxon are open
    my $is_start_open = ($fluxon->{fc_start}->{label} == -1);
    my $is_end_open = ($fluxon->{fc_end}->{label} == -2);
    my $transonic_velocity = 0;

    # If the fluxon is labeled as a plasmoid, or if both ends are open, return undefined
    return undef if($fluxon->{plasmoid} || ($is_start_open + $is_end_open != 1));


    (my $r_vr_scaled, my $r_fr_scaled, my $theta, my $phi) =
            gen_fluxon_schonflow_phys($fluxon, $distance_array_degrees, $fluxon_id);

    # Return the final radial velocity and flux expansion
    return ($r_vr_scaled, $r_fr_scaled, $theta, $phi);
}

sub gen_fluxon_schonflow_phys {
    my $me = shift;
    my $distance_array_degrees = shift;
    my $fluxon_id = shift;
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
        our $th = acos($z/$r1*$opt{r0});
        our $ph = atan2($y, $x);
        our $A = pdl(map { $_->{A}} ($me->vertices))->(-1:0:-1) * $opt{r0} * $opt{r0};
        our $bfield = pdl(map {$_->{b_vec}} ($me->vertices))->(-1:0:-1);
        our $bmag = pdl(map {$_->{b_mag}} ($me->vertices))->(-1:0:-1);
        our $T = pdl(map {$_->{T}} ($me->vertices))->(-1:0:-1);
        our $rho = pdl(map {$_->{rho}} ($me->vertices))->(-1:0:-1);

        $A(0:-2) .= $A(1:-1);
    }

    # Get rid of end anomalies
    our $A;
    $A->((0)) .= $A->((1));
    $A->((-1)) .= $A->((-2));

    # Calculate magnetic field magnitude
    our $bfield;
    our $bmag;

    # print "bfield = $bfield\n\n";
    # print "Bmag = $bmag\n\n";

    our $r1;
    our $r;
    our $th;
    our $ph;
    our $T;
    our $rho;


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

    sub wind_speed {
        my ($fss, $theta_b) = @_;

        # Define constants (From Schonfeld 2022)
        my $c1 = 2/9;
        my $c2 = 0.8;
        my $c3 = 2.0;
        my $c4 = 2.0;
        my $c5 = 3.0;

        my $v0 = 285; #km/s
        my $vm = 910 - $v0; #km/s

        my $speed = $v0 + ($vm / (1 + $fss)**$c1) * (1 - $c2 * exp(-($theta_b / $c3)**$c4))**$c5;
        return $speed;

    }

    # Define an expansion factor
    # dist_interp
    # our $fr1 = $A / $A->((1)) * ($r->((1))) * ($r->((1))) / $r / $r;
    my $Arat = $A / $A->((1));
    my $Rrat = ($r->((1))) * ($r->((1))) / $r / $r;
    our $fr1 = $Arat * $Rrat;
    # our $fr1 = pdl(($r0*$r0)/($r*$r) * $bmag->((1))/$bmag);

    # # Replace NaNs with zeros
    # $fr->inplace->setvaltobad(0)->setnantobad;
    # $fr->inplace->badmask(0);
    # # Check if all elements in $fr are zeros
    # if (all($fr == 0)) {
    #     $fr = $fr1;
    # }

    # my $th0_ind = int($th->((0)) / 3.14159265358 * 180) + 90;
    # my $ph0_ind = int($ph->((0)) / 3.14159265358 * 180) + 180;

    my $phi0 = $ph->((1));
    my $sin_theta0 = sin($th->((1)));
    my $sin_theta = sin($th);

    # print "th = " . $th->((0)) . ", \nph = " . $ph->((0)) . "\n";
    # print "sin_theta = $sin_theta\n\n";
    # print $distance_array_degrees . "\n";
    my $distance_degrees = interpolate_2d_lonlat($distance_array_degrees, $phi0, $sin_theta0);
    # print "Interpolated value: " .  . "\n\n\n";



    # my $distance = $dist_interp->(pdl($ph), pdl($sin_theta));


    # die;

    # my $ch_map = rfits($ch_map_path);
    # my ($theta_bound, $phi_bound, $distance) = gradient_descent($ch_map, $th0_ind, $ph0_ind);

    # $distance = $distance * 3.14159265358 / 180;
    # plot_image_with_points($ch_map_path, [$th0_ind, $ph0_ind], [$theta_bound, $phi_bound], $fluxon_id);


    # my $flow_field;
    # my $distance = 3;
    # print($fr1);
    my $flow_field = wind_speed($fr1, $distance_degrees/3);
    # my $flow_field = random($flow_field->dim(0)) * 500 + 250;
    print "\nDistance = $distance_degrees\n";
    print "r = $r\n";
    print "expansion factor = $fr1\n";
    print "Flow Field = $flow_field\n\n";
    # die;
    # Create a Gaussian distribution around the mean
    my $mean = 550;
    my $std_dev = 75;
    # srand(int($CR)+int($fluxon_id));
    # $flow_field = (grandom($rho->dim(0)) * $std_dev) + $mean;

    ## Calculate Return Values ##

    # This array is (r, v) in units of (r_sun, km/s)
    my $r_v_scaled = pdl($r / $opt{r0}, $flow_field);

    # This array is (r, fr) in units of (r_sun, unitless)
    my $r_fr_scaled = pdl($r / $opt{r0}, $fr1);


    # Print the PDLs and their lengths
    # print "r_v_scaled = $r_v_scaled\n";
    # print "\n\nLength of r_v_scaled: " . $r_v_scaled->nelem . "\n";
    # print "Length of r_fr_scaled: " . $r_fr_scaled->nelem . "\n\n";
    # print "r_fr_scaled = $r_fr_scaled\n";

    # Return the constructed arrays
    return ($r_v_scaled, $r_fr_scaled, $th, $ph);
}


use PDL;
use PDL::NiceSlice;

sub gradient_descent {
    my ($image, $latitude, $longitude) = @_;
    my $threshold = 0.1;

    # Convert geographic coordinates to image coordinates
    my $row = ($latitude);   # Assuming latitude and longitude are already in pixel coordinates
    my $col = ($longitude);

    # Check initial value
    my $initial_value = $image->at($row, $col);

    # Determine whether to perform descent or ascent
    my $is_ascent = $initial_value < $threshold;

    # Gradient descent or ascent
    my $max_steps = 500;  # Maximum number of steps to prevent infinite loops
    my $step_size = 1;    # Initial step size in pixels

    for (my $step = 0; $step < $max_steps; $step++) {
        my ($best_value, $best_row, $best_col) = ($initial_value, $row, $col);

        foreach my $r_offset (-$step_size, 0, $step_size) {
            foreach my $c_offset (-$step_size, 0, $step_size) {
                next if $r_offset == 0 && $c_offset == 0;  # Skip the current point

                my $new_row = $row + $r_offset;
                my $new_col = $col + $c_offset;

                # Ensure new coordinates are within bounds
                next if $new_row < 0 || $new_row >= $image->dim(0);
                next if $new_col < 0 || $new_col >= $image->dim(1);

                my $new_value = $image->at($new_row, $new_col);

                # Update best value found depending on descent or ascent
                if (($is_ascent && $new_value > $best_value) || (!$is_ascent && $new_value < $best_value)) {
                    ($best_value, $best_row, $best_col) = ($new_value, $new_row, $new_col);
                }
            }
        }

        # # Dynamically adjust the step size
        # if (abs($best_value - $initial_value) < 0.01) {
        #     $step_size += 1;
        # } else {
        #     $step_size = 1;  # Reset step size to initial value
        # }

        # Check if the best value found crosses the threshold
        if ((!$is_ascent && $best_value < $threshold) || ($is_ascent && $best_value >= $threshold)) {
            my $distance = sqrt(($best_row - $row)**2 + ($best_col - $col)**2);
            return ($best_row, $best_col, $distance);
        }

        # Update row and col for the next iteration
        $initial_value = $best_value;
        ($row, $col) = ($best_row, $best_col);
    }

    # Return NaN if threshold not crossed within max steps
    return ('NaN', 'NaN', 'NaN');
}


use PDL;
use PDL::Graphics::Gnuplot;

sub plot_image_with_points {
    my ($image_path, $original_point, $discovered_point, $fluxon_id) = @_;
    our $fh;
    # print $fluxon_id;
    # Save data and image path to a CSV file
    if ($fluxon_id == 0){
        open $fh, '>', 'data.csv' or die "Cannot open data.csv: $!";
        print $fh "x1,y1,x2,y2,image_path=$image_path\n";
    } else {
        open $fh, '>>', 'data.csv' or die "Cannot open data.csv: $!";
    }
    print $fh "$original_point->[1],$original_point->[0],$discovered_point->[1],$discovered_point->[0]\n";
    close $fh;


    # print "Data saved to data.csv\n";


}

sub do_image_plot{
    # Call the Python script for plotting
    my $python_script = 'plot_distance.py';
    my $exit_status = system("python $python_script");
    if ($exit_status == 0) {
        print "Python script executed successfully.\n";
    } else {
        print "Python script failed to execute.\n";
    }

    # <STDIN>;
}
# our $count = 0;

1;
__END__

