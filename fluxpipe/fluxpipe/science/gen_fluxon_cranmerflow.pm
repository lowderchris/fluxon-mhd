=head1 NAME

gen_fluxon_wsaflow - Comprehensive Module for Solar Wind Solution and Fluxon Flow Analysis

=cut

package gen_fluxon_wsaflow;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(gen_fluxon_cranmerflow);

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



sub gen_fluxon_cranmerflow {
    my $fluxon = shift;

    # Check if the start and end    points of the fluxon are open
    my $is_start_open = ($fluxon->{fc_start}->{label} == -1);
    my $is_end_open = ($fluxon->{fc_end}->{label} == -2);

    # If the fluxon is labeled as a plasmoid, or if both ends are open, return undefined
    return undef if($fluxon->{plasmoid} || ($is_start_open + $is_end_open != 1));


    (my $r_vr_scaled, my $r_fr_scaled, my $theta, my $phi) =
            gen_fluxon_cranmerflow_phys($fluxon, $distance_array_degrees, $fluxon_id);

    # Return the final fluxon array, fluxon radius, and magnetic field components
    return ($r_vr_scaled, $r_fr_scaled, $theta, $phi);
}

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

sub gen_fluxon_cranmerflow_phys {
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
        our $th = acos($z/$r1*$opt{r0});
        our $ph = atan2($y, $x);
        our $A = pdl(map { $_->{A}} ($me->vertices))->(-1:0:-1);
        our $bfield = pdl(map {$_->{b_vec}} ($me->vertices))->(-1:0:-1);
        our $bmag = pdl(map {$_->{b_mag}} ($me->vertices))->(-1:0:-1);
        our $T = pdl(map {$_->{T}} ($me->vertices))->(-1:0:-1);
        our $rho = pdl(map {$_->{rho}} ($me->vertices))->(-1:0:-1);

    }

    # Get rid of end anomalies
    our $A;
    # $A(0:-2) .= $A(1:-1);
    # $A->((0)) .= $A->((1));
    # $A->((-1)) .= $A->((-2));

    # Calculate magnetic field magnitude
    our $bfield;
    our $bmag;

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

        if (1){
            # Define constants (From Schonfeld 2022)
            my $c1 = 2/9;
            my $c2 = 0.8;
            my $c3 = 4.0; # 2
            my $c4 = 2.0;
            my $c5 = 3.0;

            my $v0 = 285; #km/s
            my $vm = 910 - $v0; #km/s
            our $speed = $v0 + ($vm / (1 + $fss)**$c1) * (1 - $c2 * exp(-($theta_b / $c3)**$c4))**$c5;

        } else {

            # Define constants (From Wiengarten 2014)
            my $c1 = 2/9;
            my $c2 = 0.8;
            my $c3 = 3.0; # degrees
            my $c4 = 2.0;
            my $c5 = 3.0;

            my $v0 = 200; #km/s
            my $vm = 675; #km/s
            our $speed = $v0 + ($vm / (1 + $fss)**$c1) * (1 - $c2 * exp(-($theta_b / $c3)**$c4))**$c5;
        }
        our $speed;
        return pdl($speed);

    }


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
    my $rn = $r1/$opt{r0};


    # print("The height is $height" . "\n");

    # die;

    # Find the index of the value closest to 2.5 Rs
    my $target = 2.5;
    my $differences = abs($rn - $target);
    my ($ss_ind) = $differences->minimum_ind;

    # print("Difference = $differences\n");
    # print("ss_ind = $ss_ind\n");

    # Find the expansion factor at that index
    # my $rss = interpolate_1d($rn, $rn, $target);


    my $lowind = 2;
    my $r_s = $rn->(($lowind));
    my $b_s = $bmag->(($lowind));
    my $a_s = $A->(($lowind));
    my $An = $A/$a_s;

    # my $rss = $r1->(($ss_ind));
    # my $bss = $bmag->(($ss_ind));
    # my $ass = $A->(($ss_ind));

    my $rss = $target;
    my $bss = interpolate_1d($rn, $bmag, $target);
    my $ass = interpolate_1d($rn, $An, $target);
    my $fss =   abs($b_s / $bss) * ($r_s * $r_s) / ($rss * $rss);
    my $fr_all= abs($b_s / $bmag)* ($r_s * $r_s) /  ($rn * $rn);

    my $len = $fr_all->nelem;
    # $fr_all->index(0) .= $fr_all->index(2);
    # $fr_all->index(1) .= $fr_all->index(2);
    # $fr_all->index($len-1) .= $fr_all->index($len-2);
    # $fr_all->index($len-2) .= $fr_all->index($len-3);

    # Find the distance from the coronal hole boundary at that index

    # my $phi0 = $ph->(($lowind));
    $len = $ph->nelem - 1;
    my $phi0 = pdl($ph->at($lowind));
    my $phi1 = pdl($ph->at($len));
    my $r00 = pdl($rn->at($lowind));
    my $r11 = pdl($rn->at($len));
    my $sin_theta0 = pdl(sin($th)->at($lowind));
    my $sin_theta1 = pdl(sin($th)->at($len));
    if ($phi0 < 0) {$phi0 += 2 * 3.1415926;}
    if ($phi1 < 0) {$phi1 += 2 * 3.1415926;}

    my $distance_degrees = interpolate_2d_lonlat($distance_array_degrees, $phi0, $sin_theta0);

    # Calculate the wind speed at the given location
    # my $flow_field = wind_speed($An, $distance_degrees);
    # my $speed = interpolate_1d($rn, $flow_field, 5.0);

    my $speed = wind_speed($ass, $distance_degrees);

    # write the $fid, $fss, and $distance_degrees to a file
    open my $fh, '>>', 'seq_ass_theta.csv' or die "Cannot open data.csv: $!";
    print $fh "$fluxon_id, $ass, $distance_degrees, $speed\n";
    close $fh;

    # print "\nfss = $fss\n";
    # print "Distance = $distance_degrees\n";
    # print "Flow Field = $flow_field\n\n";

    use strict;
    use warnings;
    use PDL;
    use PDL::Graphics::Gnuplot;

    # Sample data
    my $ones = ones($r->dims);

    # # Plot the data with labels
    # my $plot = PDL::Graphics::Gnuplot->new(persist => 1);
    # $plot->plot({title=>$fluxon_id}, {legend=>"phi"},$rn, $ph);
    # $plot->replot({legend=>"phi", with=>"points"},$rn, $ph);
    # $plot->replot(legend=>"ph0", with=>"points", pointsize=>3, $r00, $phi0);
    # $plot->replot(legend=>"ph0", with=>"points", pointsize=>3, $r11, $phi1);

    # $plot->replot({legend=>"stheta"},$rn, sin($th));
    # $plot->replot({legend=>"stheta", with=>"points"},$rn, sin($th));
    # $plot->replot(legend=>"sth0", with=>"points", pointsize=>3 ,$r00, $sin_theta0);
    # $plot->replot(legend=>"sth0", with=>"points", pointsize=>3 ,$r11, $sin_theta1);
    # $plot->replot({legend=>"distance"},$rn, $ones * $distance_degrees);

    # my $plot = PDL::Graphics::Gnuplot->new(persist => 1);
    # my $top_speed = $flow_field->((-1));
    # $plot->plot({title=>"ID: $fluxon_id, Top Speed: $top_speed", logscale => 'xy'}, {legend=>"rn", with=>"points"},$rn, $rn);
    # $plot->replot({legend=>"bmag", with =>'points'},$rn, $bmag);
    # $plot->replot({legend=>"fss"},$rn, $ones * $fss);
    # $plot->replot({legend=>"fr_all"},$rn, $fr_all);
    # $plot->replot({legend=>"An"},$rn, $An);
    # $plot->replot({legend=>"A"},$rn, $A);
    # $plot->replot({legend=>"flow"},$rn, $flow_field);
    # $plot->replot({legend=>"flow", with=>"points"},$rn, $flow_field);
    # $plot->replot({legend=>"flow", with=>"points"},$rn->((-1)), $flow_field->((-1)));


    # Block execution
    # print "Press ENTER to continue...\n";
    # <STDIN>;


    ## Calculate Return Values ##
    my $speed_tall = $ones * $speed;
    # This array is (r, v) in units of (r_sun, km/s)
                                            #To consider plotting:
                                                # $phi0, $sin_theta0, $distance_degrees, $flow_field
    # our $r_v_scaled = pdl($rn, $ones * ($flow_field->at($len)));
    our $r_v_scaled = pdl($rn, $speed_tall);
    # if (1-$en_open) {
    #     # print $r_v_scaled;
    # } else {
    #     our $r_v_scaled = pdl($rn, $ones * $flow_field->at(0));
    #     # print $r_v_scaled;
    # }

    # our $r_v_scaled;
    # This array is (r, fr) in units of (r_sun, unitless)
    my $r_fr_scaled = pdl($rn, $A);
    # print $fr_all;

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


    # # # Define a normalized expansion factor
    # my $Arat = $A / $A->((1));
    # my $Rrat = ($r->((1))) * ($r->((1))) / $r / $r;
    # my $Brat = $bmag / $bmag->((1));
    # our $fr1 = $Arat * $Rrat;

    # # Initialize Gnuplot object
    # use PDL::Graphics::Gnuplot;

    # gpwin('pngcairo')->plot(
    #     { title => 'Field vs. r', xlabel => 'r', ylabel => 'Field Value', logscale => 'xy' },
    #     with => 'lines', $r1/$opt{r0}, $Arat,
    #     with => 'lines', $r1/$opt{r0}, $fr1,
    #     # with => 'lines', $r, $r,
    #     with => 'lines', $r1/$opt{r0}, $bmag,
    # );

    # my $ch_map = rfits($ch_map_path);
    # my ($theta_bound, $phi_bound, $distance) = gradient_descent($ch_map, $th0_ind, $ph0_ind);

    # $distance = $distance * 3.14159265358 / 180;
    # plot_image_with_points($ch_map_path, [$th0_ind, $ph0_ind], [$theta_bound, $phi_bound], $fluxon_id);